import logging.config
import os.path


LOGGING_CONFIG = os.path.join(os.path.dirname(__file__), "logging.ini")
logging.config.fileConfig(LOGGING_CONFIG)

logger = logging.getLogger()
