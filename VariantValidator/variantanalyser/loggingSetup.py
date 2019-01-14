
# Import python diagnostic tools
import logging
from StringIO import StringIO
import traceback


# Set up logging
def loggingSetup():
    print("THIS CODE IS CALLED FUCKING ONCE")
    if "VV" in logging.Logger.manager.loggerDict:
        return

    VALIDATOR_DEBUG = os.environ.get('VALIDATOR_DEBUG')
    if VALIDATOR_DEBUG is None:
        VALIDATOR_DEBUG="info console" #Set default value
    logger=logging.getLogger("VV")
    #Set logging urgency levels.
    if "debug" in VALIDATOR_DEBUG:
        logLevel =logging.DEBUG
    elif "warning" in VALIDATOR_DEBUG:
        logLevel =logging.WARNING
    elif "info" in VALIDATOR_DEBUG:
        logLevel =logging.INFO
    elif "error" in VALIDATOR_DEBUG:
        logLevel =logging.ERROR
    elif "critical" in VALIDATOR_DEBUG:
        logLevel =logging.CRITICAL

    if "file" in VALIDATOR_DEBUG:
        logFileHandler=logging.FileHandler("VV-log.txt")
        logFileHandler.setLevel(logLevel)
        logger.addHandler(logFileHandler)
    if "console" in VALIDATOR_DEBUG:
        logConsoleHandler=logging.StreamHandler()
        logConsoleHandler.setLevel(logLevel)
        logger.addHandler(logConsoleHandler)
    #Create a log string to add to validations.
    logString=StringIO()
    logStringHandler=logging.StreamHandler(logString)
    #We want the validation metadata to not contain debug info which may change with program operation
    logStringHandler.setLevel(logging.INFO)
    logger.addHandler(logStringHandler)
    logger.setLevel(logging.DEBUG) #The logger itself must be set with an appropriate level of urgency.

    #print(logger.handers)
    logger.propagate=False

loggingSetup()

#Test
#logger.debug("Message D")
#logger.info("Message I")
#logger.warning("Message W")
#logger.error("Message E")
#logger.critical("Message C")#

#print("TEST "+logString.getvalue())

