from os import path

from luna.util.default_values import DOCK6_PATH

import logging
logger = logging.getLogger()


class Params(dict):

    def __init__(self, params=None):
        params = params or {}
        super().__init__(params)

    @property
    def names(self):
        """list: The list of parameters."""
        return sorted(self.keys())

    def save(self, output_file):
        # Don't create the output file if the params dictionary is empty.
        if not self:
            logger.info("There are no parameters to be saved.")
            return

        with open(output_file, "w") as OUT:
            for k, v in self.items():
                k = k.ljust(61)
                OUT.write("%s%s\n" % (k, str(v)))

    @classmethod
    def load(cls, input_file):
        params = {}

        with open(input_file) as IN:
            for line in IN:
                line = line.rstrip("\n")

                name, val = line.split(" ", 1)
                val = val.lstrip(" ")

                try:
                    val = int(val)
                except ValueError:
                    try:
                        val = float(val)
                    except ValueError:
                        pass
                params[name] = val

        return cls(params)


class DefaultParams(Params):

    def __init__(self, **kwargs):
        dock_path = path.abspath(path.join(path.realpath(__file__), '../'))
        default_file_path = f"{dock_path}/default.cfg"

        params = Params.load(default_file_path)

        for k, v in kwargs.items():
            params[k] = v

        super().__init__(params)
