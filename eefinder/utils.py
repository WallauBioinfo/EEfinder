from datetime import datetime
import re
from pathlib import Path

def check_outdir(outdir: str) -> str:
    if  outdir.endswith("/"):
        outdir = re.sub("/$", "", outdir)
    
    Path(outdir).mkdir(parents=True, exist_ok=True)

    return outdir


def step_info(step: str, start_time: str, end_time: str, message: str) -> dict:
    total_time_minutes = (end_time - start_time) / 60
    start_time_formated = datetime.fromtimestamp(start_time).strftime(
        "%Y-%m-%d %H:%M:%S"
    )
    end_time_formated = datetime.fromtimestamp(end_time).strftime("%Y-%m-%d %H:%M%:%S")

    return {
        "step": step,
        "start_time": start_time_formated,
        "end_time": end_time_formated,
        "total_time_minutes": f"{total_time_minutes:.4f}",
        "message": message,
    }


def running_info(
    arguments: list, start_time: str, end_time: str, steps_infos: dict
) -> dict:
    total_time_minutes = (end_time - start_time) / 60
    start_time_formated = datetime.fromtimestamp(start_time).strftime(
        "%Y-%m-%d %H:%M:%S"
    )
    end_time_formated = datetime.fromtimestamp(end_time).strftime("%Y-%m-%d %H:%M%:%S")
    arguments_info = {
        "genome_file": arguments[0],
        "prefix": arguments[15],
        "outdir": arguments[1],
        "database": arguments[2],
        "dbmetadata": arguments[3],
        "baits": arguments[4],
        "mode": arguments[5],
        "length": arguments[6],
        "flank": arguments[7],
        "limit": arguments[8],
        "range_junction": arguments[9],
        "mask_per": arguments[10],
        "clean_masked": arguments[11],
        "threads": arguments[12],
        "removetmp": arguments[13],
        "index_databases": arguments[14],
        "merge_level": arguments[6],
    }

    return {
        "arguments": arguments_info,
        "start_time": start_time_formated,
        "end_time": end_time_formated,
        "total_time_minutes": f"{total_time_minutes:.4f}",
        "steps_information": steps_infos,
    }
