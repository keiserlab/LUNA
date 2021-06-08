import configparser
import logging

logger = logging.getLogger()


class Config:
    def __init__(self, conf_file):
        self.config = configparser.ConfigParser()
        try:
            self.config.read(conf_file)
        except Exception as e:
            logger.exception(e)
            raise IOError('Configuration file %s not read.' % conf_file)

    def get_section_map(self, section):
        section_map = {}
        options = self.config.options(section)
        for option in options:
            try:
                section_map[option] = self.config.get(section, option)
            except Exception as e:
                logger.exception(e)
                section_map[option] = None
        return section_map

    def __getattr__(self, attr):
        if hasattr(self.config, attr):
            return getattr(self.config, attr)
        else:
            raise AttributeError("The attribute '%s' does not exist in the class %s." % (attr, self.__class__.__name__))
