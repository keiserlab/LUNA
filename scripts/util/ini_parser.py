import configparser
import logging

logger = logging.getLogger()


class Config():

    def __init__(self, iniFile):
        self.config = configparser.ConfigParser()
        try:
            self.config.read(iniFile)
        except Exception as e:
            logger.exception(e)
            raise IOError("Configuration file %s not read." % iniFile)

    def get_section_map(self, section):
        sectionMap = {}
        options = self.config.options(section)
        for option in options:
            try:
                sectionMap[option] = self.config.get(section, option)
            except Exception as e:
                logger.exception(e)
                logger.warning("Exception on option '%s'!" % option)
                sectionMap[option] = None
        return sectionMap
