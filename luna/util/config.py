import configparser
import ast
import logging

logger = logging.getLogger()


class Config:
    """
    Parser for configuration files.

    Parameters
    ----------
    config_file : str
        The pathname for the configuration file.

    Attributes
    ----------

    config : `configparser.ConfigParser`
        The parsed configuration file.
    """

    def __init__(self, config_file):
        self.config = configparser.ConfigParser(allow_no_value=True)
        try:
            self.config.read(config_file)
        except Exception as e:
            logger.exception(e)
            raise IOError('Configuration file %s not read.' % config_file)

    def parse_value(self, section, property_name):
        """Transform the value of property ``property_name`` from
        section ``section`` to a Python object.

        Returns
        -------
         : object
        """
        value = self.config.get(section, property_name)
        return ast.literal_eval(value)

    def get_section_map(self, section):
        """Try to access the section ``section`` from the parsed
        configuration file.

        Returns
        -------
         : dict
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

    def as_dict(self, only_props=True):
        """Return the object parameters as :py:class:`dict`."""
        if only_props:
            return {k: v
                    for s in self.config.sections()
                    for k, v in self.config.items(s)}
        return {s: dict(self.config.items(s)) for s in self.config.sections()}

    def __getattr__(self, attr):
        if hasattr(self.config, attr):
            return getattr(self.config, attr)
        else:
            raise AttributeError("The attribute '%s' does not exist in the "
                                 "class %s." % (attr,
                                                self.__class__.__name__))
