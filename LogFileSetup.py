import logging

def configure_logging(logfile_path):
    '''
    Initialize logging defaults for in-notebook messages and 'log.info' file written to store intermediate results during analysis of each defect
    '''

    # Set default format for each line of log messages within notebook
    notebook_formatter = logging.Formatter("[%(levelname)s] [Cell line num: %(lineno)s] %(message)s")
    # Set default format for each line in log.info file (look into methods to outputt cell num, not just line num in cell)
 #   info_file_formatter = logging.Formatter("[%(levelname)s] [Notebook cell num: %(???)s] [Cell line num: %(lineno)s] %(message)s")

    # Initialise log.info for defect processing information
    defect_logger = logging.getLogger()
    # For log.info file
    info_file_handler = logging.FileHandler(logfile_path + ".info", mode='w')
    info_file_handler.setLevel(logging.INFO)
#    info_file_handler.setFormatter(info_file_formatter)
    # For messages within notebook
    notebook_handler = logging.StreamHandler()
    notebook_handler.setLevel(logging.INFO)
    notebook_handler.setFormatter(notebook_formatter)

    # Remove default handlers and add custom ones (for log.info file and messages in notebooks)
    list(map(defect_logger.removeHandler, defect_logger.handlers[:]))
    list(map(defect_logger.removeFilter, defect_logger.filters[:]))
    defect_logger.setLevel(logging.INFO)
    defect_logger.addHandler(info_file_handler)
    defect_logger.addHandler(notebook_handler)

    return defect_logger
