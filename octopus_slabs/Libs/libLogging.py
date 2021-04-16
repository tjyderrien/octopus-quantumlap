# https://gist.githubusercontent.com/agalea91/8be90216a184e24f54123f5922d635b7/raw/841d7bc7c4806dd00722c9d04085d916873f9b81/init_logger.py

import logging
import colorlog

def init_logger(dunder_name, verbose) -> logging.Logger:

    #         '%(asctime)s - '
    log_format = (
        '%(name)s - '
        '%(funcName)s - '
        '%(levelname)s - '
        '%(message)s'
    )
    bold_seq = '\033[1m'
    colorlog_format = (
        f'{bold_seq} '
        '%(log_color)s '
        f'{log_format}'
    )
    colorlog.basicConfig(format=colorlog_format)
    logger = logging.getLogger(dunder_name)

    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # Output full log
    fh = logging.FileHandler('app.log')
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter(log_format)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Output warning log
    fh = logging.FileHandler('app.warning.log')
    fh.setLevel(logging.WARNING)
    formatter = logging.Formatter(log_format)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Output error log
    fh = logging.FileHandler('app.error.log')
    fh.setLevel(logging.ERROR)
    formatter = logging.Formatter(log_format)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger
