import logging


def new_logging_file(filename):
    # Set the LOG file at the working path
    fh = logging.FileHandler(filename, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - '
                                  '%(levelname)s - %(filename)s - '
                                  '%(funcName)s - line %(lineno)d => '
                                  '%(message)s')
    fh.setFormatter(formatter)

    logger = logging.getLogger(__name__)

    # Remove the existing file handlers
    for hdlr in logger.handlers[:]:
        if isinstance(hdlr, logger.FileHander):
            logger.removeHandler(hdlr)
    # Set the new handler
    logger.addHandler(fh)
    # Set the log level to INFO, DEBUG as the default is ERROR
    logger.setLevel("WARNING")

    return logger
