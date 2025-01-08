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
    """Create and configure a logger with specified level and file path.

    Sets up a global logger with a StreamHandler and basic configuration. The
    logger is named "autogrow" and will output to both console and file.

    Args:
        level (int): Logging level (e.g., logging.DEBUG, logging.INFO)
        file_path (str): Path where log file will be created. Defaults to
            "log.txt"

    Note:
        This function modifies global logger, handler, and log_filename
        variables.
    """
    global logger
    global handler
    global log_filename

    log_filename = file_path

    # Create a logger
    logger = logging.getLogger("autogrow")
    logger.setLevel(level)

    # Create a handler
    handler = logging.StreamHandler()

    # Add the handler to the logger
    # logger.addHandler(handler)
    logging.basicConfig(level=logging.DEBUG, handlers=[handler])

    set_log_tab_level(0)


class CustomFormatter(logging.Formatter):
    """Custom logging formatter that uses abbreviated level names.

    Customizes the format of log messages by using shorter level names (e.g.,
    "DEBG" instead of "DEBUG") while maintaining timestamp and message
    formatting.

    Args:
        fmt (str, optional): Message format string. Defaults to timestamp-level-
            message format
        datefmt (str, optional): Date format string. Defaults to ISO-like format
        style (str, optional): Style of format string. Defaults to "%"
    """

    def __init__(self, fmt=None, datefmt=None, style="%"):
        """
        Initialize the formatter with the specified format and date format.
        
        Args:
            fmt (str, optional): Message format string. Defaults to
                timestamp-level-message format
            datefmt (str, optional): Date format string. Defaults to ISO-like
                format
            style (str, optional): Style of format string. Defaults to "%"
        """
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
        """
        Format the log record with the custom level name.
        
        Args:
            record (logging.LogRecord): The log record to format
        
        Returns:
            str: The formatted log message
        """
        # Replace the levelname with our custom string
        record.levelname = self.level_strings.get(record.levelno, record.levelname)
        return super().format(record)


def set_log_tab_level(tab_count: int):
    """
    Set the indentation level for log messages.

    Configures a new formatter with the specified indentation level and applies
    it to the global handler.

    Args:
        tab_count (int): Number of tab indentations to prepend to log messages

    Note:
        Modifies the global handler's formatter and level variables.
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
    """Wrap a message to fit within the maximum line length.

    Args:
        msg (str): The message to wrap

    Returns:
        List[str]: List of strings, each representing a wrapped line of the
            original message

    Note:
        Uses MAX_MSG_LINE_LENGTH (60) as the maximum width for wrapping.
    """
    return textwrap.wrap(msg, width=MAX_MSG_LINE_LENGTH, break_long_words=True)


def increase_log_tab_level():
    """Increases the indentation level for log messages by one tab.

    Increments the global indentation level and updates the log formatter
    accordingly.

    Note:
        Modifies the global level variable.
    """
    global level
    set_log_tab_level(level + 1)


def decrease_log_tab_level():
    """Decreases the indentation level for log messages by one tab.

    Decrements the global indentation level if it's greater than zero and
    updates the log formatter accordingly.

    Note:
        Modifies the global level variable. Will not decrease level below zero.
    """
    global level
    if level > 0:
        set_log_tab_level(level - 1)


def log_info(msg: str):
    """Log a message at INFO level with proper indentation and line wrapping.

    Writes the message to both the logger (if configured) and the log file. Long
    messages are wrapped to fit MAX_MSG_LINE_LENGTH, with subsequent lines
    indented.

    Args:
        msg (str): The message to log

    Note:
        Uses global level and logger variables. Writes to global log_filename.
        Wrapped lines are indented with INDENT.
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
    """Log a message at DEBUG level with proper indentation and line wrapping.

    Writes the message to both the logger (if configured) and the log file. Long
    messages are wrapped to fit MAX_MSG_LINE_LENGTH, with subsequent lines
    indented.

    Args:
        msg (str): The message to log

    Note:
        Uses global level and logger variables. Writes to global log_filename.
        Wrapped lines are indented with INDENT.
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
    Log a message at WARNING level with proper indentation and line wrapping.

    Writes the message to both the logger (if configured) and the log file. Long
    messages are wrapped to fit MAX_MSG_LINE_LENGTH, with subsequent lines
    indented.

    Args:
        msg (str): The message to log

    Note:
        Uses global level and logger variables. Writes to global log_filename.
        Wrapped lines are indented with INDENT.
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
    """Context manager for temporarily adjusting log indentation levels.

    Increases the log indentation level upon entering the context and decreases
    it upon exit, ensuring proper nesting of logged messages.

    Example:
        >>> with LogLevel():
        ...     # Messages logged here will be indented one level
        ...     log_info("This message is indented")
    """

    def __enter__(self):
        """Increases log indentation level when entering the context.

        Returns:
            LogLevel: The context manager instance
        """
        increase_log_tab_level()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Decreases log indentation level when exiting the context.

        Args:
            exc_type: Type of any exception that occurred
            exc_value: Instance of any exception that occurred
            traceback: Traceback of any exception that occurred
        """
        decrease_log_tab_level()


class WithoutLogging:
    """Context manager that temporarily disables all logging operations.

    Provides a context where all logging is disabled, including for newly
    created loggers. Restores the original logging state upon exit.

    Example:
        >>> with WithoutLogging():
        ...     # All logging operations in this block will be suppressed
        ...     logger.info("This won't be logged")
    """

    def __init__(self):
        """Initialize context manager by storing the original logging state."""
        self.original_logging_class = logging.getLoggerClass()
        self.original_get_logger = logging.getLogger

    def disabled_getLogger(self, name=None):
        """
        Get a disabled logger instance.

        Args:
            name (str, optional): Name for the logger. Defaults to None

        Returns:
            logging.Logger: A disabled logger instance
        """
        logger = self.original_get_logger(name)
        logger.disabled = True
        return logger

    def __enter__(self):
        """Disable all existing and future loggers when entering the context."""
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
        """Restore logging functionality when exiting the context."""
        # Re-enable all loggers
        for logger in logging.Logger.manager.loggerDict.values():
            if isinstance(logger, logging.Logger):
                logger.disabled = False
        logging.root.disabled = False

        # Restore original logger class and getLogger function
        logging.setLoggerClass(self.original_logging_class)
        logging.getLogger = self.original_get_logger
