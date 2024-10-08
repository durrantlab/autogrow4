"""
Module to set up and manage logging.

It provides functions to create loggers, set logging levels, and manage
indentation levels for structured logging output.
"""

import logging
import os
from typing import Optional

logger: Optional[logging.Logger] = None
handler: Optional[logging.StreamHandler] = None
log_filename = "log.txt"
level = 0


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
    logger.addHandler(handler)

    set_log_tab_level(0)


def set_log_tab_level(tab_count: int):
    """
    Set the log tab level, which affects the indentation of log messages.

    Args:
        tab_count (int): The number of tabs to use.
    """
    global level

    # Create an initial formatter
    new_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - " + (tab_count * "  ") + "%(message)s"
    )

    if handler is not None:
        handler.setFormatter(new_formatter)

    level = tab_count


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

    if logger is not None:
        # Add space to align with DEBUG
        logger.info(f" {msg}")

    with open(log_filename, "a") as f:
        f.write((level * "  ") + msg + "\n")


def log_debug(msg: str):
    """
    Log a debug message.

    Args:
        msg (str): The message to log.
    """
    global level
    global logger

    if logger is not None:
        logger.debug(msg)

    with open(log_filename, "a") as f:
        f.write((level * "  ") + msg + "\n")


def log_warning(msg: str):
    """
    Log a warning message.

    Args:
        msg (str): The message to log.
    """
    global level
    global logger

    if logger is not None:
        logger.warning(msg)

    with open(log_filename, "a") as f:
        f.write((level * "  ") + msg + "\n")


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
