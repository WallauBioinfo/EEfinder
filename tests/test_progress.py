"""Unit tests for the progress display and the download retry/stall logic."""

from __future__ import annotations

import io
import sys
import time

import pytest

from eefinder import progress
from eefinder.progress import (
    DEFAULT_ATTEMPTS,
    NO_PROGRESS_ENV,
    STALLED,
    progress_bar,
    run_streaming,
    run_with_retries,
    show_progress,
)


class _Tty(io.StringIO):
    """A stream that claims to be a terminal."""

    def isatty(self) -> bool:
        return True


class TestShowProgress:
    def test_hidden_when_not_a_terminal(self, monkeypatch):
        monkeypatch.delenv(NO_PROGRESS_ENV, raising=False)
        assert not show_progress(io.StringIO())

    def test_shown_on_a_terminal(self, monkeypatch):
        monkeypatch.delenv(NO_PROGRESS_ENV, raising=False)
        assert show_progress(_Tty())

    def test_the_environment_variable_wins(self, monkeypatch):
        monkeypatch.setenv(NO_PROGRESS_ENV, "1")
        assert not show_progress(_Tty())

    def test_an_empty_or_zero_value_does_not_disable_it(self, monkeypatch):
        monkeypatch.setenv(NO_PROGRESS_ENV, "0")
        assert show_progress(_Tty())


class TestProgressBar:
    def test_every_item_is_yielded_without_a_terminal(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        with progress_bar([1, 2, 3], "Working") as items:
            assert list(items) == [1, 2, 3]

    def test_every_item_is_yielded_with_a_bar(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: True)
        with progress_bar(iter([1, 2, 3]), "Working", length=3) as items:
            assert list(items) == [1, 2, 3]


class TestRunStreaming:
    def test_output_is_captured_and_the_status_returned(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        assert run_streaming(["echo", "hello"]) == (0, "hello")

    def test_a_failure_keeps_the_tools_own_message(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        returncode, output = run_streaming(["sh", "-c", "echo boom >&2; exit 3"])
        assert (returncode, output) == (3, "boom")

    def test_output_is_mirrored_when_progress_is_shown(self, monkeypatch, capsys):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: True)
        run_streaming(["echo", "visible"])
        assert "visible" in capsys.readouterr().err

    def test_only_the_tail_is_kept(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        _, output = run_streaming(["sh", "-c", "printf 'abcdefghij'"], keep_chars=4)
        assert output == "ghij"


class TestStallDetection:
    def test_a_silent_process_is_killed(self, monkeypatch):
        """NCBI transfers hang with the connection open; no status ever reports it."""
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        start = time.monotonic()
        returncode, _ = run_streaming(["sleep", "30"], stall_timeout=1)
        assert returncode == STALLED
        assert time.monotonic() - start < 10

    def test_a_growing_output_file_counts_as_life(self, monkeypatch, tmp_path):
        """A quiet but working transfer must not be mistaken for a hang."""
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        target = tmp_path / "download.zip"
        returncode, _ = run_streaming(
            [
                "sh",
                "-c",
                f"for i in 1 2 3 4; do printf x >> {target}; sleep 0.6; done",
            ],
            stall_timeout=1.5,
            liveness_path=str(target),
        )
        assert returncode == 0

    def test_no_timeout_means_wait(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        assert run_streaming(["sleep", "0.2"], stall_timeout=None)[0] == 0


class TestRunWithRetries:
    def test_a_transient_failure_is_retried(self, monkeypatch, tmp_path):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        monkeypatch.setattr(progress, "RETRY_BACKOFF", 0)
        marker = tmp_path / "attempts"
        # Fails the first time, succeeds once the marker exists.
        command = [
            "sh",
            "-c",
            f"if [ -f {marker} ]; then exit 0; else touch {marker}; exit 1; fi",
        ]
        assert run_with_retries(command, attempts=3)[0] == 0

    def test_the_last_failure_is_returned(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        monkeypatch.setattr(progress, "RETRY_BACKOFF", 0)
        returncode, output = run_with_retries(
            ["sh", "-c", "echo nope >&2; exit 7"], attempts=2
        )
        assert (returncode, output) == (7, "nope")

    def test_the_partial_output_is_discarded_between_attempts(
        self, monkeypatch, tmp_path
    ):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        monkeypatch.setattr(progress, "RETRY_BACKOFF", 0)
        cleaned = []
        run_with_retries(
            ["sh", "-c", "exit 1"],
            attempts=3,
            before_retry=lambda attempt: cleaned.append(attempt),
        )
        assert cleaned == [2, 3]

    def test_a_single_attempt_disables_retrying(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        calls = []

        def fake_streaming(command, **kwargs):
            calls.append(command)
            return 1, ""

        monkeypatch.setattr(progress, "run_streaming", fake_streaming)
        run_with_retries(["true"], attempts=1)
        assert len(calls) == 1

    def test_stalled_attempts_are_retried_too(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        monkeypatch.setattr(progress, "RETRY_BACKOFF", 0)
        outcomes = [(STALLED, "stalled"), (0, "fine")]

        monkeypatch.setattr(
            progress, "run_streaming", lambda command, **kwargs: outcomes.pop(0)
        )
        assert run_with_retries(["true"], attempts=3)[0] == 0
        assert not outcomes


def test_default_attempts_is_three():
    """The documented default: NCBI fails intermittently, three tries is enough."""
    assert DEFAULT_ATTEMPTS == 3


class TestPermanentErrors:
    """A wrong request fails identically every time; only transfers are retried."""

    @pytest.mark.parametrize(
        "message",
        [
            "Error: The taxonomy name 'NotARealTaxon' is not a recognized virus taxon.",
            "The taxonomy ID '0' does not match any existing taxids",
            "unknown flag: --nope",
        ],
    )
    def test_request_errors_are_permanent(self, message):
        assert progress.is_permanent_error(message)

    @pytest.mark.parametrize(
        "message",
        [
            "Error: Internal error (invalid zip archive). Please try again",
            "connection reset by peer",
            "giving up after 1 attempt(s)",
            "",
        ],
    )
    def test_transfer_errors_are_worth_retrying(self, message):
        assert not progress.is_permanent_error(message)

    def test_a_permanent_failure_is_not_retried(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        monkeypatch.setattr(progress, "RETRY_BACKOFF", 0)
        calls = []

        def fake_streaming(command, **kwargs):
            calls.append(command)
            return 1, "Error: The taxonomy name 'Nope' is not a recognized virus taxon."

        monkeypatch.setattr(progress, "run_streaming", fake_streaming)
        run_with_retries(["true"], attempts=3)
        assert len(calls) == 1


class _TimedStream(io.StringIO):
    """A terminal-like stream that records when each write happened."""

    def __init__(self):
        super().__init__()
        self.writes = []

    def isatty(self) -> bool:
        return True

    def write(self, text):
        if text:
            self.writes.append((time.monotonic(), text))
        return super().write(text)


class TestOutputIsForwardedImmediately:
    """Regression: a partial line must not wait for a full buffer to flush.

    A buffered text read of a fixed size blocks until that many characters
    exist, so a tool drawing a progress bar appeared frozen mid-word until
    enough further output arrived.
    """

    def test_a_short_message_appears_while_the_process_is_still_running(
        self, monkeypatch
    ):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: True)
        stream = _TimedStream()
        monkeypatch.setattr(sys, "stderr", stream)

        start = time.monotonic()
        # Six characters, far below any read buffer, then a long pause.
        run_streaming(["sh", "-c", "printf Valida; sleep 2"], stall_timeout=None)
        finished = time.monotonic()

        assert stream.writes, "nothing was mirrored"
        first_write, text = stream.writes[0]
        assert "Valida" in text
        # It reached the terminal about when it was printed, not at exit.
        assert first_write - start < 1.0
        assert finished - first_write > 1.0

    def test_multibyte_output_is_not_corrupted(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        _, output = run_streaming(["sh", "-c", "printf 'progresso: 100%% ✓'"])
        assert output.endswith("✓")


class TestExitIsDetectedWithoutEndOfFile:
    """Regression: a lingering grandchild must not keep the run hanging.

    A command can exit while a background process it started still holds the
    write end of the pipe open. Waiting for end-of-file then hangs *after* the
    tool has printed everything and finished -- which is what a stuck
    `Validating package files ... 100%` looked like.
    """

    def test_the_run_ends_when_the_command_exits(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        start = time.monotonic()
        # The shell exits immediately; the backgrounded sleep inherits stdout and
        # keeps the pipe open for 10s.
        returncode, output = run_streaming(
            ["sh", "-c", "sleep 10 & printf finished"], stall_timeout=None
        )
        elapsed = time.monotonic() - start

        assert returncode == 0
        assert "finished" in output
        assert elapsed < 5, f"waited {elapsed:.1f}s for a command that had exited"

    def test_output_written_just_before_exit_is_still_captured(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        _, output = run_streaming(
            ["sh", "-c", "sleep 5 & printf 'Validating 100%%'; exit 0"],
            stall_timeout=None,
        )
        assert "Validating 100%" in output


class TestStandardInputIsClosed:
    """Regression: an inherited stdin can keep a finished command alive forever.

    The NCBI `datasets` client completes its download, prints its validation bar
    and then never exits when it inherits an open stdin -- the symptom being a
    run apparently frozen on `Validating package files ... 100%`. Measured on the
    same download: 3 seconds with stdin closed, indefinite with it open.
    """

    @pytest.mark.parametrize("mirrored", [False, True])
    def test_a_command_reading_stdin_still_finishes(self, monkeypatch, mirrored):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: mirrored)
        start = time.monotonic()
        # `cat` reads stdin until end-of-file; with stdin left open it would
        # never return.
        returncode, output = run_streaming(
            ["sh", "-c", "printf 'Validating 100%%'; cat"], stall_timeout=None
        )
        elapsed = time.monotonic() - start

        assert returncode == 0
        assert "Validating 100%" in output
        assert elapsed < 5, f"the command waited {elapsed:.1f}s on stdin"


class TestRetryEnvironment:
    """A failure with a known workaround is retried *differently*."""

    def test_the_environment_returned_is_applied_to_later_attempts(
        self, monkeypatch, tmp_path
    ):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        monkeypatch.setattr(progress, "RETRY_BACKOFF", 0)
        marker = tmp_path / "attempt"
        # Fails until GODEBUG is set, which is what the fallback provides.
        command = [
            "sh",
            "-c",
            'if [ -n "$GODEBUG" ]; then printf ok; exit 0; fi; '
            "printf 'stream error: INTERNAL_ERROR' >&2; exit 1",
        ]
        returncode, output = run_with_retries(
            command,
            attempts=3,
            retry_env=lambda out: (
                {"GODEBUG": "http2client=0"} if "stream error" in out else None
            ),
        )
        assert (returncode, output) == (0, "ok")
        assert not marker.exists()

    def test_no_environment_is_applied_when_the_hook_declines(self, monkeypatch):
        monkeypatch.setattr(progress, "show_progress", lambda stream=None: False)
        monkeypatch.setattr(progress, "RETRY_BACKOFF", 0)
        seen = []
        run_with_retries(
            ["sh", "-c", 'printf "%s" "${GODEBUG:-unset}" >&2; exit 1'],
            attempts=2,
            retry_env=lambda out: seen.append(out) or None,
        )
        assert seen  # the hook was consulted
