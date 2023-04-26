import logging
import logging.config
import os
import colorlog
import gzip


FILE_FORMAT = '[%(asctime)s]    %(levelname)-8s %(filename)16s:%(lineno)-10d %(processName)-20s %(message)s'
CONSOLE_FORMAT = '[%(asctime)s]    %(log_color)s%(levelname)-10s %(reset)s%(filename)16s:%(lineno)-10d %(message)s'

VERBOSITY_LEVEL = {4: logging.DEBUG,
                   3: logging.INFO,
                   2: logging.WARNING,
                   1: logging.ERROR,
                   0: logging.CRITICAL}


class CompressedRotatingFileHandler(logging.handlers.RotatingFileHandler):

    # Source: http://roadtodistributed.blogspot.com/2011/04/compressed-rotatingfilehandler-for.html
    def doRollover(self):
        self.stream.close()
        if self.backupCount > 0:
            for i in range(self.backupCount - 1, 0, -1):
                sfn = "%s.%d.gz" % (self.baseFilename, i)
                dfn = "%s.%d.gz" % (self.baseFilename, i + 1)

                if os.path.exists(sfn):
                    if os.path.exists(dfn):
                        os.remove(dfn)
                    os.rename(sfn, dfn)
            dfn = self.baseFilename + ".1.gz"

            if os.path.exists(dfn):
                os.remove(dfn)
            try:
                f_in = open(self.baseFilename, 'rb')
                f_out = gzip.open(dfn, 'wb')
                f_out.writelines(f_in)
            except Exception:
                if not os.path.exists(dfn):
                    if os.path.exists(self.baseFilename):
                        os.rename(self.baseFilename, dfn)
            finally:
                if "f_out" in dir() and f_out is not None:
                    f_out.close()
                if "f_in" in dir() and f_in is not None:
                    f_in.close()
            if os.path.exists(self.baseFilename):
                os.remove(self.baseFilename)
        self.mode = 'w'
        self.stream = self._open()


def new_logging_file(filename, logging_level=logging.INFO,
                     logger_name=None, propagate=True,
                     log_format=None, mode='a'):
    # Set the LOG file at the working path
    fh = CompressedRotatingFileHandler(filename, mode, 5 * 1024 * 1024, 20)
    file_format = log_format or FILE_FORMAT
    file_formatter = logging.Formatter(file_format, '%Y-%m-%d %H:%M:%S')
    fh.setFormatter(file_formatter)
    fh.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    console_format = log_format or CONSOLE_FORMAT
    console_formatter = \
        colorlog.ColoredFormatter(console_format, '%Y-%m-%d %H:%M:%S',
                                  log_colors={'DEBUG': 'cyan',
                                              'INFO': 'green',
                                              'WARNING': 'yellow',
                                              'ERROR': 'red',
                                              'CRITICAL': 'bold_red'})
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


def load_default_logging_config():
    logger = logging.getLogger()

    # Remove the existing file handlers
    for hdlr in logger.handlers[:]:
        logger.removeHandler(hdlr)

    LOGGING_CONFIG = os.path.join(os.path.dirname(__file__), "logging.cfg")
    logging.config.fileConfig(LOGGING_CONFIG)

    return logger
