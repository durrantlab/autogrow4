"""
Module to set up and manage logging.

It provides functions to create loggers, set logging levels, and manage
indentation levels for structured logging output.
"""
import logging
import os
from typing import List, Optional
import textwrap

logger: Optional[logging.Logger] = None
handler: Optional[logging.StreamHandler] = None
log_filename = "log.txt"
level = 0
MAX_MSG_LINE_LENGTH = 60
TAB = "    "
INDENT = " "


def create_logger(level: int, file_path: str = "log.txt"):
    """
    Create and configure the logger with the specified level.

    Args:
        level (int): Logging level (e.g., logging.DEBUG, logging.INFO).
        file_path (str): Dated results directory.
    """
    global logger
    global handler
    global log_filename

    log_filename = file_path

    # Create a logger
    logger = logging.getLogger("grid_ml")
    logger.setLevel(level)

    # Create a handler
    handler = logging.StreamHandler()

    # Add the handler to the logger
    # logger.addHandler(handler)
    logging.basicConfig(level=logging.DEBUG, handlers=[handler])

    set_log_tab_level(0)


class CustomFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style="%"):
        if fmt is None:
            fmt = "%(asctime)s - %(levelname)s - %(message)s"
        if datefmt is None:
            datefmt = "%Y-%m-%d %H:%M:%S"
        super().__init__(fmt, datefmt, style)
        self.level_strings = {
            logging.DEBUG: "DEBG",
            logging.INFO: "INFO",
            logging.WARNING: "WARN",
            logging.ERROR: "ERRO",
            logging.CRITICAL: "CRIT",
        }

    def format(self, record):
        # Replace the levelname with our custom string
        record.levelname = self.level_strings.get(record.levelno, record.levelname)
        return super().format(record)


def set_log_tab_level(tab_count: int):
    """
    Set the log tab level, which affects the indentation of log messages.

    Args:
        tab_count (int): The number of tabs to use.
    """
    global level, handler

    # Create an initial formatter
    new_formatter = CustomFormatter(
        fmt="%(asctime)s - %(levelname)s - " + (tab_count * TAB) + "%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if handler is not None:
        handler.setFormatter(new_formatter)

    level = tab_count


def wrap_msg(msg: str) -> List[str]:
    """
    Wrap a message to fit within the maximum line length using textwrap.

    Args:
        msg (str): The message to wrap.

    Returns:
        List[str]: The wrapped message.
    """
    return textwrap.wrap(msg, width=MAX_MSG_LINE_LENGTH, break_long_words=True)


def increase_log_tab_level():
    """
    Increase the log tab level, which affects the indentation of log messages.
    """
    global level
    set_log_tab_level(level + 1)


def decrease_log_tab_level():
    """
    Decrease the log tab level, which affects the indentation of log messages.
    Ensures that the level does not go below zero.
    """
    global level
    if level > 0:
        set_log_tab_level(level - 1)


def log_info(msg: str):
    """
    Log an informational message.

    Args:
        msg (str): The message to log.
    """
    global level
    global logger

    for i, m in enumerate(wrap_msg(msg)):
        if i > 0:
            m = INDENT + m

        if logger is not None:
            # Add space to align with DEBUG
            logger.info(m)

        with open(log_filename, "a") as f:
            f.write((level * TAB) + m + "\n")


def log_debug(msg: str):
    """
    Log a debug message.

    Args:
        msg (str): The message to log.
    """
    global level
    global logger

    for i, m in enumerate(wrap_msg(msg)):
        if i > 0:
            m = INDENT + m
        if logger is not None:
            logger.debug(m)

        with open(log_filename, "a") as f:
            f.write((level * TAB) + m + "\n")


def log_warning(msg: str):
    """
    Log a warning message.

    Args:
        msg (str): The message to log.
    """
    global level
    global logger

    for i, m in enumerate(wrap_msg(msg)):
        if i > 0:
            m = INDENT + m
        if logger is not None:
            logger.warning(m)

        with open(log_filename, "a") as f:
            f.write((level * TAB) + m + "\n")


class LogLevel:
    """Context manager to increase and decrease the log tab level."""

    def __enter__(self):
        """
        Increase the log tab level when entering a block of code.
        """
        increase_log_tab_level()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Decrease the log tab level when exiting a block of code.
        """
        decrease_log_tab_level()


class WithoutLogging:
    """Context manager to temporarily disable all logging, including new loggers."""

    def __init__(self):
        self.original_logging_class = logging.getLoggerClass()
        self.original_get_logger = logging.getLogger

    def disabled_getLogger(self, name=None):
        """Custom getLogger function that returns disabled loggers."""
        logger = self.original_get_logger(name)
        logger.disabled = True
        return logger

    def __enter__(self):
        """
        Disable all existing and future loggers when entering the context.
        """
        # Disable all existing loggers
        for logger in logging.Logger.manager.loggerDict.values():
            if isinstance(logger, logging.Logger):
                logger.disabled = True
        logging.root.disabled = True

        # Create a new Logger subclass that's always disabled
        class DisabledLogger(self.original_logging_class):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.disabled = True

        # Override logger creation to ensure new loggers are disabled
        logging.setLoggerClass(DisabledLogger)
        logging.getLogger = self.disabled_getLogger

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Restore logging functionality when exiting the context.
        """
        # Re-enable all loggers
        for logger in logging.Logger.manager.loggerDict.values():
            if isinstance(logger, logging.Logger):
                logger.disabled = False
        logging.root.disabled = False

        # Restore original logger class and getLogger function
        logging.setLoggerClass(self.original_logging_class)
        logging.getLogger = self.original_get_logger
