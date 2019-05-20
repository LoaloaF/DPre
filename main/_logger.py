"""Definition of the logger"""
import logging

# default logger
logger = logging.getLogger('dpre')
stream_handler = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s-DPre %(funcName)s:\n%(message)s')
stream_handler.setFormatter(formatter)
logger.setLevel(logging.INFO)
logger.addHandler(stream_handler)

# spacer is used to make free lines
spacer = logging.getLogger('DPre')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(logging.Formatter(''))
spacer.setLevel(logging.INFO)
spacer.addHandler(stream_handler)

log_init = '=|=|=|=|=|=|=|=|====INITIATION====|=|=|=|=|=|=|=|=|='
log_plot = '=|=|=|=|=|=|=|=|====PLOT====|=|=|=|=|=|=|=|=|='