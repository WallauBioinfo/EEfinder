import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("eefinder")


def enable_debug() -> None:
    """Lower the ``eefinder`` logger to DEBUG.

    Called by the ``--debug`` CLI flag of ``screening`` / ``get-databases`` so
    that ``logger.debug(...)`` messages are emitted; without it the logger stays
    at the INFO level configured above and debug messages are suppressed.
    """
    logger.setLevel(logging.DEBUG)
