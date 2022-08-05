from os import path

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

    def __init__(self, config=None):
        config = config or {}
        super().__init__(config)

    @property
    def params(self):
        """list: The list of parameters."""
        return sorted(self.keys())

    @classmethod
    def from_config_file(cls, config_file):

        if not path.exists(config_file):
            raise OSError("File '%s' does not exist." % config_file)

        config = {}
        params = Config(config_file)
        for section in params.sections():
            params_dict = params.get_section_map(section)
            params_dict = {prop: params.parse_value(section, prop)
                           for prop in params_dict}
            config.update(params_dict)

        return cls(config)


class DefaultInteractionConfig(InteractionConfig):
    """Default parameters for calculating interactions in LUNA."""

    def __init__(self, **kwargs):

        file_path = path.abspath(path.join(path.realpath(__file__), '../'))
        default_config_file = f"{file_path}/config.cfg"

        config = InteractionConfig.from_config_file(default_config_file)

        for k, v in kwargs.items():
            config[k] = v

        super().__init__(config)
