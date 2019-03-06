import logging, logging.config

def configure_logging(logfile_path):
    '''
    Initialize logging defaults for 'log.info' file written to store metadata during analysis of each defect
    '''

    default_formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] [%(name)s] [%(funcName)s():%(lineno)s] %(message)s")
    root_logger = logging.getLogger()
    info_file_handler = logging.FileHandler(logfile_path + ".info", mode='w')
    info_file_handler.setLevel(logging.INFO)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(default_formatter)
    list(map(root_logger.removeHandler, root_logger.handlers[:]))
    list(map(root_logger.removeFilter, root_logger.filters[:]))
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(info_file_handler)
    root_logger.addHandler(console_handler)

    return root_logger
