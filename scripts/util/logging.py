import logging


DEFAULT_FORMAT = '%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(filename)s - %(funcName)s - line %(lineno)d => %(message)s'


def new_logging_file(filename, log_format=DEFAULT_FORMAT, mode='a'):
    # Set the LOG file at the working path
    fh = logging.FileHandler(filename, mode)
    formatter = logging.Formatter(log_format)
    fh.setFormatter(formatter)

    logger = logging.getLogger()

    # Remove the existing file handlers
    for hdlr in logger.handlers[:]:
        if isinstance(hdlr, logging.FileHandler):
            logger.removeHandler(hdlr)
    # Set the new handler
    logger.addHandler(fh)
    # Set the log level to INFO, DEBUG as the default is ERROR
    logger.setLevel("INFO")
