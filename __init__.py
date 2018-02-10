import logging

outfilePath='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/'

log_file = outfilePath+'runlog.txt'
log_level = logging.DEBUG
logging.basicConfig(filename=log_file, level=log_level, format='%(asctime)s %(message)s')
