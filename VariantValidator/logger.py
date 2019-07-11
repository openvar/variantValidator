import logging.config
from configparser import ConfigParser
from . import settings

# Change settings based on config
config = ConfigParser()
config.read(settings.CONFIG_DIR)

if config['logging'].getboolean('log') is True:
    settings.LOGGING_CONFIG['handlers']['console']['level'] = config['logging']['console'].upper()
    settings.LOGGING_CONFIG['handlers']['file']['level'] = config['logging']['file'].upper()

    logging.config.dictConfig(settings.LOGGING_CONFIG)
else:
    logging.getLogger('VariantValidator').addHandler(logging.NullHandler())

