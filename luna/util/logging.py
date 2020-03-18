import logging
import logging.config
import os.path


DEFAULT_FORMAT = '%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(filename)s - %(funcName)s - line %(lineno)d => %(message)s'


def new_logging_file(filename, logger_name=None, propagate=True, log_format=DEFAULT_FORMAT, mode='a'):
    # Set the LOG file at the working path
    fh = logging.handlers.RotatingFileHandler(filename, mode, 5 * 1024 * 1024, 100)
    formatter = logging.Formatter(log_format)
    fh.setFormatter(formatter)
    fh.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    ch.setLevel(logging.ERROR)

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
