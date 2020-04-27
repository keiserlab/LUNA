import logging
import logging.config
import os.path
import colorlog


FILE_FORMAT = '[%(asctime)s]    %(levelname)-8s %(filename)16s:%(lineno)-10d %(threadName)-16s %(message)s'
CONSOLE_FORMAT = '[%(asctime)s]    %(log_color)s%(levelname)-10s %(reset)s%(filename)16s:%(lineno)-10d %(message)s'


def new_logging_file(filename, logging_level=logging.INFO, logger_name=None, propagate=True, log_format=None, mode='a'):
    # Set the LOG file at the working path
    fh = logging.handlers.RotatingFileHandler(filename, mode, 5 * 1024 * 1024, 100)
    file_format = log_format or FILE_FORMAT
    file_formatter = logging.Formatter(file_format, '%Y-%m-%d %H:%M:%S')
    fh.setFormatter(file_formatter)
    fh.setLevel(logging_level)

    ch = logging.StreamHandler()
    console_format = log_format or CONSOLE_FORMAT
    console_formatter = colorlog.ColoredFormatter(console_format, '%Y-%m-%d %H:%M:%S',
                                                  log_colors={'DEBUG': 'cyan',
                                                              'INFO': 'green',
                                                              'WARNING': 'yellow',
                                                              'ERROR': 'red',
                                                              'CRITICAL': 'red,bg_white'})
    ch.setFormatter(console_formatter)
    ch.setLevel(logging_level)

    logger = logging.getLogger(logger_name)

    if logger_name is not None:
        logger.propagate = propagate

    # Remove the existing file handlers
    for hdlr in logger.handlers[:]:
        logger.removeHandler(hdlr)

    # Set the new handlers
    logger.addHandler(fh)
    logger.addHandler(ch)


def load_default_logging_conf():
    logger = logging.getLogger()

    # Remove the existing file handlers
    for hdlr in logger.handlers[:]:
        logger.removeHandler(hdlr)

    LOGGING_CONF = os.path.join(os.path.dirname(__file__), "logging.ini")
    logging.config.fileConfig(LOGGING_CONF)
