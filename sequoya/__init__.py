import logging

logger = logging.getLogger('Sequoya')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

# deactivate some loggers
logging.getLogger('pyMSA').setLevel(logging.CRITICAL)
logging.getLogger('jMetalPy').setLevel(logging.CRITICAL)
