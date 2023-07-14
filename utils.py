from datetime import datetime
from pathlib import Path
import json
import logging
import numpy as np


logger = logging.getLogger()

_ticks = []
_prefixes = []


def tic(message: str = ""):
    _ticks.append(datetime.now())
    _prefixes.append(message)
    return


def toc(message: str = "", pop: bool = True):
    now = datetime.now()
    tick = _ticks.pop() if pop else _ticks[-1]
    prefix = _prefixes.pop() if pop else _prefixes[-1]
    logger.info(f"{now-tick} - {' '*4*len(_ticks)}{prefix}{message}")
    return


def timeit(label=None):
    def inner(func):
        def wrapper(*args, **kwargs):
            tic(f"{label if label else func.__name__}()")
            ret = func(*args, **kwargs)
            toc(" (total)")
            return ret

        return wrapper

    return inner


def read_json(file):
    with open(file) as f:
        json_data = json.load(f)
    return json_data


# From https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.float32):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


# TODO: Compare with new pandas function to_json()
def save_json(
    filename: str, json_data: object, indent: int = None, separators=None
) -> None:
    with open(filename, "w") as handle:
        json.dump(
            json_data, handle, cls=NumpyEncoder, indent=indent, separators=separators
        )
    return


def setup_logging(level: str = "INFO", directory: Path = None):
    levels = {
        "DEBUG": logging.DEBUG,  # 10
        "INFO": logging.INFO,  # 20
        "WARN": logging.WARN,  # 30
        "WARNING": logging.WARNING,  # 30
        "ERROR": logging.ERROR,  # 40
        "CRITICAL": logging.CRITICAL,  # 50
        "FATAL": logging.FATAL,  # 50
    }

    # logger = logging.getLogger()    # did this at module scope above
    logger.setLevel(levels[level])
    console = logging.StreamHandler()
    console.setLevel(levels[level])
    formatter = logging.Formatter("%(levelname)s:%(message)s")
    console.setFormatter(formatter)
    logger.addHandler(console)
    disk = logging.FileHandler(directory / "logfile.txt")
    disk.setLevel(levels[level])
    disk.setFormatter(formatter)
    logger.addHandler(disk)

    return
