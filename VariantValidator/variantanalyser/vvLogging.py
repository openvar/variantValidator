
import logging
import os
from StringIO import StringIO


class logger():
    #Grand unified variant validator logging static class.
    logger=logging.getLogger("VV")
    logString=""
    @staticmethod
    def loggingSetup():
        # Set up logging
        # I need to use the VVObfuscator in the logger global dictionary
        # becuase it's a global variable tied to the logger module
        # Modules are singletons, but their variables are not. Consequently
        # this is the only sensible way to ensure that the logging setup is called
        # once. If another programmer has any better ideas that leave these functions
        # with a configured VV logger object that only has its handlers added once,
        # feel free to fix it up.
        if "VVObfuscator" in logging.Logger.manager.loggerDict:
            return
        logging.getLogger("VVObfuscator")
        # Check envrionment variables
        VALIDATOR_DEBUG = os.environ.get('VALIDATOR_DEBUG')
        if VALIDATOR_DEBUG is None:
            VALIDATOR_DEBUG = "info console"  # Set default value
        logger.logger = logging.getLogger("VV")
        # Set logging urgency levels.
        if "debug" in VALIDATOR_DEBUG:
            logLevel = logging.DEBUG
        elif "warning" in VALIDATOR_DEBUG:
            logLevel = logging.WARNING
        elif "info" in VALIDATOR_DEBUG:
            logLevel = logging.INFO
        elif "error" in VALIDATOR_DEBUG:
            logLevel = logging.ERROR
        elif "critical" in VALIDATOR_DEBUG:
            logLevel = logging.CRITICAL

        if "file" in VALIDATOR_DEBUG:
            logFileHandler = logging.FileHandler("VV-log.txt")
            logFileHandler.setLevel(logLevel)
            logger.logger.addHandler(logFileHandler)
        if "console" in VALIDATOR_DEBUG:
            logConsoleHandler = logging.StreamHandler()
            logConsoleHandler.setLevel(logLevel)
            logger.logger.addHandler(logConsoleHandler)
        # Create a log string to add to validations.
        logger.logString = StringIO()
        logStringHandler = logging.StreamHandler(logger.logString)
        # We want the validation metadata to not contain debug info which may change with program operation
        logStringHandler.setLevel(logging.INFO)
        logger.logger.addHandler(logStringHandler)
        logger.logger.setLevel(logging.DEBUG)  # The logger itself must be set with an appropriate level of urgency.

        # print(logger.handers)
        logger.logger.propagate = False
    @staticmethod
    def debug(s):
        logger.logger.debug("DEBUG: "+s)
    @staticmethod
    def info(s):
        logger.logger.info("INFO : "+s)
    @staticmethod
    def warning(s):
        logger.logger.warn("WARN : "+s)
    @staticmethod
    def error(s):
        logger.logger.error("ERROR: "+s)
    @staticmethod
    def critical(s):
        logger.logger.critical("CRIT : "+s)
    @staticmethod
    def trace(s):
        logger.logger.debug("TRACE: "+s)
    @staticmethod
    def getString(s):
        return logger.logString.getvalue()


logger.loggingSetup()



#Test
#logger.debug("Message D")
#logger.info("Message I")
#logger.warning("Message W")
#logger.error("Message E")
#logger.critical("Message C")#

#print("TEST "+logString.getvalue())

