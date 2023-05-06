from os import path
from collections import defaultdict

from luna.util.config import Config

import logging
logger = logging.getLogger()


class InteractionConfig(dict):
    """Generic class to define parameters for interactions.

    Parameters
    ----------
    config : dict, optional
        A dict containing parameters for calculating interactions.

    Examples
    --------

    >>> config = {"max_ha_dist_hb_inter": 2.5}
    >>> inter_config = InteractionConfig(config)
    >>> print(inter_config["max_ha_dist_hb_inter"])
    2.5
    inter_config["max_ha_dist_hb_inter"] = 3
    >>> print(inter_config["max_ha_dist_hb_inter"])
    3
    """

    def __init__(self, config=None, names_map=None):
        config = config or {}
        self.names_map = names_map or {}
        super().__init__(config)

    @property
    def params(self):
        """list: The list of parameters."""
        return sorted(self.keys())

    @classmethod
    def from_config_file(cls, config_file):
        """Initialize from a configuration file.

        Parameters
        ----------
        config_file : str
            The configuration file pathname.

        Returns
        -------
         : `InteractionConfig`
        """
        if not path.exists(config_file):
            raise OSError("File '%s' does not exist." % config_file)

        config = {}
        params = Config(config_file)
        names_map = {}

        for section in params.sections():
            params_dict = params.get_section_map(section)
            params_dict = {prop: params.parse_value(section, prop)
                           for prop in params_dict}
            names_map.update({k: section for k in params_dict})
            config.update(params_dict)

        return cls(config, names_map)

    def save_config_file(self, config_file):
        """Save the interaction parameters into a configuration file.

        Parameters
        ----------
        config_file : str
            The output configuration file.
        """
        args_by_section = defaultdict(list)
        for k, v in self.items():
            section = self.names_map.get(k, "Other")
            args_by_section[section].append((k, v))

        with open(config_file, "w") as OUT:
            for section in sorted(args_by_section):
                OUT.write("[%s]\n" % section)
                for arg, val in sorted(args_by_section[section]):
                    OUT.write("%s = %s\n" % (arg, val))
                OUT.write("\n")


class DefaultInteractionConfig(InteractionConfig):
    """Default parameters for calculating interactions in LUNA."""

    def __init__(self, **kwargs):

        file_path = path.abspath(path.join(path.realpath(__file__), '../'))
        default_config_file = f"{file_path}/config.cfg"

        config = InteractionConfig.from_config_file(default_config_file)

        for k, v in kwargs.items():
            config[k] = v

        super().__init__(config)
