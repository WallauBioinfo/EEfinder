"""Progress reporting for the long-running steps of ``get-databases``.

Downloading and preparing a RefSeq database takes minutes to hours, and until
now EEfinder printed nothing while it happened: the ``datasets`` CLI writes a
perfectly good progress display, but it was being swallowed by
``capture_output=True``. This module fixes that and adds bars for the steps
EEfinder itself performs.

Two rules keep the output honest:

* progress is shown **only on a terminal**. Piped into a file or a workflow
  engine, the escape sequences that redraw a progress bar are noise, so the
  behaviour falls back to what it was before (silent, output captured);
* mirroring a subprocess's output does not stop EEfinder from reading it -- the
  tail is still kept so a failure can be reported with the tool's own message.

Set ``EEFINDER_NO_PROGRESS=1`` to turn the display off even on a terminal.

A command is judged by **what it produced**, not by whether it exited: when a
``success_check`` reports the expected artifact complete, the command is stopped
and the run continues. The NCBI ``datasets`` client writes the archive, prints
its validation bar and then keeps running; waiting for it froze the pipeline on a
download that had actually finished.

Every command is run with its **standard input closed**. This is not a detail:
the NCBI ``datasets`` client finishes its download, prints its completed
validation bar and then *never exits* if it inherits an open stdin -- which is
exactly what happens in an interactive terminal. Measured on the same download,
same machine: 3 seconds with stdin closed, indefinite with it open. Nothing here
feeds a command through stdin, so closing it is free.

The same machinery detects a **stalled** download. NCBI transfers sometimes hang
with the connection open, which no exit status ever reports: the process simply
never finishes. :func:`run_with_retries` watches two liveness signals -- new
output from the tool, and growth of the file it is writing -- and when both stay
quiet for ``stall_timeout`` seconds it kills the attempt and starts over.
"""

from __future__ import annotations

import codecs
import os
import queue
import subprocess
import sys
import threading
import time

import click

from eefinder.log import logger

#: Environment variable that suppresses every progress display.
NO_PROGRESS_ENV = "EEFINDER_NO_PROGRESS"

#: How much of a mirrored subprocess's output is kept for error reporting.
KEEP_CHARS = 4000

#: Seconds of complete silence -- no new output, no growth of the output file --
#: after which a download is treated as stalled rather than slow.
DEFAULT_STALL_TIMEOUT = 180

#: How many times a download is attempted before giving up.
DEFAULT_ATTEMPTS = 3

#: Seconds to wait before a retry (multiplied by the attempt number).
RETRY_BACKOFF = 5

#: Exit status used to report a stalled attempt (no real exit status exists).
STALLED = -9

#: Messages that mean retrying is pointless: the request itself is wrong, so the
#: next attempt would fail identically. Everything else (dropped connections,
#: server errors, hangs) is worth another try.
PERMANENT_ERRORS = (
    "is not a recognized",
    "does not match any existing",
    "unknown flag",
    "unknown command",
    "no such file or directory",
)

#: Width every bar label is padded to, so the bars line up under each other.
LABEL_WIDTH = 18

#: How many times a bar is redrawn over the whole iteration.
BAR_STEPS = 200

#: Upper bound of one read from a child's pipe. Unbuffered binary reads return
#: whatever is available rather than waiting for this many bytes.
CHUNK_BYTES = 4096

#: How often the loop wakes up to check whether the command has exited.
POLL_INTERVAL = 0.2

#: Grace period for collecting output still in flight once a command exits.
DRAIN_SECONDS = 1.0

#: How long the output artifact must stay unchanged before it is examined for
#: completeness. Guards against inspecting a file that is still being written.
ARTIFACT_SETTLE = 2.0

#: How often to say that a quiet command is still being waited on. Silence is
#: indistinguishable from a freeze otherwise, and the user cannot tell whether
#: anything will ever happen.
QUIET_NOTICE_INTERVAL = 20


def show_progress(stream=None) -> bool:
    """Whether progress should be displayed on ``stream`` (default: stderr)."""
    if os.environ.get(NO_PROGRESS_ENV, "").strip() not in ("", "0", "false", "False"):
        return False
    stream = stream or sys.stderr
    try:
        return bool(stream.isatty())
    except (AttributeError, ValueError):  # pragma: no cover - closed/odd streams
        return False


def _reader(stream, sink: "queue.Queue") -> None:
    """Push a subprocess's output into ``sink`` as soon as any arrives.

    The pipe is read **unbuffered and in binary**, because a buffered text read
    of a fixed size blocks until that many characters exist: a tool that draws a
    progress bar would then appear frozen mid-word (``Valida``) until enough
    further output accumulated to fill the buffer.
    """
    try:
        while True:
            chunk = stream.read(CHUNK_BYTES)
            if not chunk:
                break
            sink.put(chunk)
    finally:
        sink.put(None)


def run_streaming(
    command: "list[str]",
    keep_chars: int = KEEP_CHARS,
    stall_timeout: "float | None" = None,
    liveness_path: "str | None" = None,
    success_check=None,
    env: "dict[str, str] | None" = None,
) -> "tuple[int, str]":
    """Run ``command``, mirroring its output and watching for a stall.

    The child's stdout and stderr are merged and forwarded byte by byte, so a
    tool that redraws its own progress bar (``datasets`` writes
    ``Downloading: virus.zip 65.5kB 125kB/s`` and updates it in place) displays
    exactly as it would when run directly. Output is mirrored only on a terminal;
    it is always read, so the last ``keep_chars`` characters can be quoted when
    the command fails.

    Parameters
    ----------
    command : list[str]
        Already-split command to execute.
    keep_chars : int
        Size of the retained tail.
    stall_timeout : float, optional
        Seconds without any sign of life after which the process is killed and
        :data:`STALLED` is returned. ``None`` waits indefinitely.
    liveness_path : str, optional
        A file the command is writing; its growth counts as a sign of life, so a
        quiet-but-working transfer is not mistaken for a hang.
    success_check : callable, optional
        Returns ``True`` once the command's actual goal is achieved -- for a
        download, that the archive on disk is complete and valid. When it does,
        the command is stopped and reported as successful **without waiting for
        it to exit**: the NCBI ``datasets`` client is known to keep running after
        it has written and validated the archive, and what EEfinder needs is the
        archive, not the exit status.
    env : dict[str, str], optional
        Extra environment variables for the child, merged over the current
        environment.

    Returns
    -------
    tuple[int, str]
        ``(exit status, captured tail)``.
    """
    mirror = show_progress()
    child_env = {**os.environ, **env} if env else None
    process = subprocess.Popen(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=0,  # binary and unbuffered: forward output the moment it arrives
        env=child_env,
    )
    decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
    chunks: "queue.Queue" = queue.Queue()
    reader = threading.Thread(
        target=_reader, args=(process.stdout, chunks), daemon=True
    )
    reader.start()

    tail = ""
    last_activity = time.monotonic()
    last_notice = last_activity
    last_size = _file_size(liveness_path)
    while True:
        try:
            chunk = chunks.get(timeout=POLL_INTERVAL)
        except queue.Empty:
            if process.poll() is not None:
                # The command finished. Do not keep waiting for end-of-file: a
                # background process it left behind can hold the write end of the
                # pipe open indefinitely, and the run would hang after the tool
                # had already printed everything and exited.
                logger.debug("command exited; draining any remaining output")
                tail = (tail + _drain(chunks, decoder, mirror))[-keep_chars:]
                break
            size = _file_size(liveness_path)
            now = time.monotonic()
            if size != last_size:
                last_size, last_activity, last_notice = size, now, now
            elif (
                success_check is not None
                and now - last_activity >= ARTIFACT_SETTLE
                and success_check()
            ):
                # The command produced what was asked of it. Whether it ever
                # returns is its own business.
                logger.debug(
                    "the expected output is complete; stopping the command "
                    "instead of waiting for it to exit"
                )
                tail = (tail + _drain(chunks, decoder, mirror))[-keep_chars:]
                _terminate(process)
                return 0, tail.strip()
            elif now - last_notice >= QUIET_NOTICE_INTERVAL:
                last_notice = now
                quiet = now - last_activity
                logger.info(
                    f"Still waiting on {os.path.basename(command[0])}: no output "
                    f"and no new data for {quiet:.0f}s"
                    + (
                        f" (giving up and retrying at {stall_timeout:.0f}s)"
                        if stall_timeout
                        else " (no stall timeout set)"
                    )
                )
            if stall_timeout and now - last_activity > stall_timeout:
                process.kill()
                process.wait()
                logger.warning(
                    f"No progress for {stall_timeout:.0f}s; the transfer looks "
                    "stalled and was interrupted."
                )
                return STALLED, tail.strip()
            continue
        if chunk is None:
            break
        last_activity = time.monotonic()
        # Decode incrementally: a chunk boundary can fall inside a multi-byte
        # character.
        text = decoder.decode(chunk)
        if not text:
            continue
        if mirror:
            sys.stderr.write(text)
            sys.stderr.flush()
        tail = (tail + text)[-keep_chars:]

    process.stdout.close()
    return process.wait(), tail.strip()


def _terminate(process: "subprocess.Popen") -> None:
    """Ask a process to stop, then insist."""
    process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()


def _drain(chunks: "queue.Queue", decoder, mirror: bool) -> str:
    """Collect whatever output is already queued, without waiting for more."""
    collected = ""
    deadline = time.monotonic() + DRAIN_SECONDS
    while time.monotonic() < deadline:
        try:
            chunk = chunks.get(timeout=0.05)
        except queue.Empty:
            break
        if chunk is None:
            break
        text = decoder.decode(chunk)
        if text and mirror:
            sys.stderr.write(text)
            sys.stderr.flush()
        collected += text
    return collected


def _file_size(path: "str | None") -> int:
    """Size of ``path`` in bytes, or ``-1`` when it does not exist yet."""
    if not path:
        return -1
    try:
        return os.path.getsize(path)
    except OSError:
        return -1


def is_permanent_error(output: str) -> bool:
    """Whether a failure is in the request rather than in the transfer.

    A misspelled taxon fails the same way every time, so retrying it only makes
    the user wait longer for the same message.
    """
    lowered = (output or "").lower()
    return any(pattern in lowered for pattern in PERMANENT_ERRORS)


def run_with_retries(
    command: "list[str]",
    attempts: int = DEFAULT_ATTEMPTS,
    stall_timeout: "float | None" = DEFAULT_STALL_TIMEOUT,
    liveness_path: "str | None" = None,
    before_retry=None,
    success_check=None,
    retry_env=None,
) -> "tuple[int, str]":
    """Run a command, retrying it when it fails or stalls.

    NCBI transfers fail in two ways that both need the same answer: an error the
    tool reports, and a hang it never reports at all. Both are retried, up to
    ``attempts`` times, with a partial output discarded in between --- except for
    failures that are permanent (:func:`is_permanent_error`), such as a
    misspelled taxon, where retrying only delays the same message.

    Parameters
    ----------
    command : list[str]
        Already-split command to execute.
    attempts : int
        Maximum number of tries (``1`` disables retrying).
    stall_timeout : float, optional
        Passed to :func:`run_streaming`.
    liveness_path : str, optional
        Passed to :func:`run_streaming`, and removed between attempts.
    before_retry : callable, optional
        Called with the attempt number just before each retry (e.g. to clean up
        a partial download).
    success_check : callable, optional
        Passed to :func:`run_streaming`: what the command was asked to produce,
        checked directly instead of trusting it to exit.
    retry_env : callable, optional
        Called with the failed attempt's output; may return environment
        variables to apply to every later attempt. This is how a failure mode
        that has a known workaround (an HTTP/2 stream reset, say) is retried
        differently instead of identically.

    Returns
    -------
    tuple[int, str]
        The ``(exit status, captured tail)`` of the last attempt.
    """
    attempts = max(1, int(attempts))
    returncode, output = 0, ""
    env: "dict[str, str]" = {}
    for attempt in range(1, attempts + 1):
        if attempt > 1:
            delay = RETRY_BACKOFF * (attempt - 1)
            if show_progress():
                # The failed attempt left the cursor mid-line and the next one
                # redraws over it; start the retry on a clean line so the warning
                # below stays readable.
                sys.stderr.write("\n")
                sys.stderr.flush()
            logger.warning(
                f"Attempt {attempt - 1} of {attempts} failed; retrying in "
                f"{delay}s ({'stalled' if returncode == STALLED else f'exit {returncode}'})"
            )
            if before_retry is not None:
                before_retry(attempt)
            time.sleep(delay)
        returncode, output = run_streaming(
            command,
            stall_timeout=stall_timeout,
            liveness_path=liveness_path,
            success_check=success_check,
            env=env or None,
        )
        if returncode == 0:
            if attempt > 1:
                logger.info(f"Succeeded on attempt {attempt} of {attempts}")
            return returncode, output
        if is_permanent_error(output):
            logger.debug("the request itself was rejected; not retrying")
            return returncode, output
        if retry_env is not None:
            extra = retry_env(output)
            if extra:
                env.update(extra)
    return returncode, output


def progress_bar(iterable, label: str, length: "int | None" = None):
    """Iterate ``iterable`` behind a progress bar when on a terminal.

    Falls back to the plain iterable otherwise, so nothing changes for piped or
    scripted runs.

    Parameters
    ----------
    iterable : iterable
        What to iterate.
    label : str
        Text shown next to the bar.
    length : int, optional
        Total number of items, for iterables without a length.

    Returns
    -------
    A context manager yielding an iterator (mirroring ``click.progressbar``).
    """
    if not show_progress():
        return _NullBar(iterable)
    return click.progressbar(
        iterable,
        label=label.ljust(LABEL_WIDTH),
        length=length,
        file=sys.stderr,
        show_pos=True,
        # A whole-RefSeq download iterates millions of records; redrawing the bar
        # for each one costs more than the work being measured.
        update_min_steps=max(1, (length or 0) // BAR_STEPS),
    )


class _NullBar:
    """A no-op stand-in for ``click.progressbar`` used off the terminal."""

    def __init__(self, iterable):
        self._iterable = iterable

    def __enter__(self):
        return self._iterable

    def __exit__(self, *exc_info):
        return False

    def update(self, _count: int = 1) -> None:
        """Accept the ``click.progressbar`` API without doing anything."""
