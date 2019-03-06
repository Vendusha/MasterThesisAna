import logging

def init_logging(logger, log_level = "INFO"):
    """ initializes logging """
    log = logging.getLogger(logger)  # get reference to logger
    # test if we have set up the logger before
    if not len(log.handlers):
        # perform setup by adding handler:
        formatter = logging.Formatter('%(asctime)s %(name)s(%(levelname)s): %(message)s',"%H:%M:%S")
        handler_stream = logging.StreamHandler()
        handler_stream.setFormatter(formatter)
        log.addHandler(handler_stream)
        # using this decorator, we can count the number of error messages
        class callcounted(object):
            """Decorator to determine number of calls for a method"""
            def __init__(self,method):
                self.method=method
                self.counter=0
            def __call__(self,*args,**kwargs):
                self.counter+=1
                return self.method(*args,**kwargs)
        log.error=callcounted(log.error)
        # set the logging level
        numeric_level = getattr(logging, log_level,
                                None)  # default: INFO messages and above
        log.setLevel(numeric_level)
    return log

