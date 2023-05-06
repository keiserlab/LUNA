import os

from luna.util.exceptions import LicenseNotFoundError


def check_license():

    if ("CHEMAXON_LICENSE_URL" not in os.environ
            and "CHEMAXON_LICENSE_SERVER_KEY" not in os.environ):
        raise LicenseNotFoundError("The Chemaxon license is not installed.")
