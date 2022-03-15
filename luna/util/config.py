import configparser
import logging

logger = logging.getLogger()


class Config:
    """
    Parser for configuration files.

    Parameters
    ----------
    conf_file : str
        The pathname for the configuration file.

    Attributes
    ----------

    config : `configparser.ConfigParser`
        The parsed configuration file.
    """

    def __init__(self, conf_file):
        self.config = configparser.ConfigParser(allow_no_value=True)
        try:
            self.config.read(conf_file)
        except Exception as e:
            logger.exception(e)
            raise IOError('Configuration file %s not read.' % conf_file)

    def get_section_map(self, section):
        """
        Try to access the section ``section`` from the parsed configuration file.
        """
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
