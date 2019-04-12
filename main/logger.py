import logging

logger = logging.getLogger('Dpre')
stream_handler = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s-DPre %(funcName)s:\n%(message)s')
stream_handler.setFormatter(formatter)
logger.setLevel(logging.INFO)
logger.addHandler(stream_handler)

spacer = logging.getLogger('DPre')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(logging.Formatter(''))
spacer.setLevel(logging.INFO)
spacer.addHandler(stream_handler)
