import logging


def new_logging_file(filename):
    # Set the LOG file at the working path
    fh = logging.FileHandler(filename, 'a')
    formatter = logging.Formatter('%(asctime)s - %(processName)s - %(name)s - '
                                  '%(levelname)s - %(filename)s - '
                                  '%(funcName)s - line %(lineno)d => '
                                  '%(message)s')
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
