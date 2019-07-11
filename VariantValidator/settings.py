import os

CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.variantvalidator')

LOG_FILE = os.path.join(os.path.expanduser('~'), '.vv_errorlog')

LOGGING_CONFIG = {
    'version': 1,
    'formatters': {
        'simple': {
            'class': 'logging.Formatter',
            'format': '%(levelname)s: %(message)s'
        },
        'detailed': {
            'class': 'logging.Formatter',
            'format': '%(asctime)s %(name)-5s %(funcName)-10s (line %(lineno)d) %(levelname)-8s %(message)s'
        }
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'level': 'DEBUG',
            'formatter': 'simple'
        },
        'file': {
            'class': 'logging.FileHandler',
            'level': 'ERROR',
            'filename': LOG_FILE,
            'mode': 'a',
            'formatter': 'detailed',
        },
    },
    'loggers': {
        'VariantValidator': {
            'level': 'DEBUG',
            'handlers': ['console', 'file'],
            'propagate': 'no',
        }
    }
}
