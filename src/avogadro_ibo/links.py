import platform
import subprocess
import logging

from .config import resolve_output_dir

logger = logging.getLogger(__name__)


def open_calcs_dir(cjson: dict) -> dict:
    calc_dir = resolve_output_dir().resolve()
    logger.debug(f"Opening calculations directory: {calc_dir}")
    if platform.system() == "Windows":
        subprocess.run(["explorer.exe", str(calc_dir)])
    elif platform.system() == "Darwin":
        subprocess.run(["open", str(calc_dir)])
    else:
        subprocess.run(["xdg-open", str(calc_dir)])

    return {
        "moleculeFormat": "cjson",
        "cjson": cjson,
    }
