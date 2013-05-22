import logging
from logging import DEBUG, INFO, WARNING, CRITICAL, ERROR
import sys


class NullHandler(logging.Handler):
    def emit(self, record):
        pass


def config_logger(name, format='%(message)s', datefmt=None,
                  stream=sys.stdout, level=logging.INFO,
                  filename=None, filemode='w', filelevel=None,
                  propagate=False):
    """
    Do basic configuration for the logging system. Similar to
    logging.basicConfig but the logger ``name`` is configurable and both a file
    output and a stream output can be created. Returns a logger object.

    The default behaviour is to create a StreamHandler which writes to
    sys.stdout, set a formatter using the "%(message)s" format string, and
    add the handler to the ``name`` logger.

    A number of optional keyword arguments may be specified, which can alter
    the default behaviour.

    :param name: Logger name
    :param format: handler format string
    :param datefmt: handler date/time format specifier
    :param stream: initialize the StreamHandler using ``stream``
            (None disables the stream, default=sys.stdout)
    :param level: logger level (default=INFO).
    :param filename: create FileHandler using ``filename`` (default=None)
    :param filemode: open ``filename`` with specified filemode ('w' or 'a')
    :param filelevel: logger level for file logger (default=``level``)
    :param propagate: propagate message to parent (default=False)

    :returns: logging.Logger object
    """
    # Get a logger for the specified name
    logger = logging.getLogger(name)
    logger.setLevel(level)
    fmt = logging.Formatter(format, datefmt)
    logger.propagate = propagate

    # Remove existing handlers, otherwise multiple handlers can accrue
    for hdlr in logger.handlers:
        logger.removeHandler(hdlr)

    # Add handlers. Add NullHandler if no file or stream output so that
    # modules don't emit a warning about no handler.
    if not (filename or stream):
        logger.addHandler(NullHandler())

    if filename:
        hdlr = logging.FileHandler(filename, filemode)
        if filelevel is None:
            filelevel = level
        hdlr.setLevel(filelevel)
        hdlr.setFormatter(fmt)
        logger.addHandler(hdlr)

    if stream:
        hdlr = logging.StreamHandler(stream)
        hdlr.setLevel(level)
        hdlr.setFormatter(fmt)
        logger.addHandler(hdlr)

    return logger
