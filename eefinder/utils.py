from datetime import datetime
import re
from pathlib import Path
import pandas as pd
from eefinder.log import logger

EXPECTED_METADATA_COLUMNS = [
    "Accession",
    "Species",
    "Genus",
    "Family",
    "Molecule_type",
    "Protein",
    "Host",
]


def check_metadata_columns(columns: list) -> list:
    """
    Check the columns of the metadata table parsed with -mt, which is read by column
    position on the taxonomy steps.

    Keyword arguments:
    columns: column names of the metadata table, in file order

    Raise a ValueError if any expected column is missing, warn if the expected columns
    are present in a different order or with extra columns, and return the expected
    columns in the order EEfinder needs them.
    """
    missing_columns = [
        column for column in EXPECTED_METADATA_COLUMNS if column not in columns
    ]
    if missing_columns:
        raise ValueError(
            f"the metadata file does not have the column(s): {', '.join(missing_columns)}. "
            f"The metadata file must have the columns: {', '.join(EXPECTED_METADATA_COLUMNS)}."
        )

    extra_columns = [
        column for column in columns if column not in EXPECTED_METADATA_COLUMNS
    ]
    if extra_columns:
        logger.warning(
            f"The metadata file has extra column(s): {', '.join(extra_columns)}. "
            "They will be ignored."
        )

    expected_columns_order = [
        column for column in columns if column in EXPECTED_METADATA_COLUMNS
    ]
    if expected_columns_order != EXPECTED_METADATA_COLUMNS:
        logger.warning(
            f"The metadata file columns are not in the expected order: {', '.join(EXPECTED_METADATA_COLUMNS)}. "
            "They will be reordered in memory, the metadata file is not modified."
        )

    return EXPECTED_METADATA_COLUMNS


def check_metadata_file(metadata_file: str) -> list:
    """
    Check the header of the metadata table parsed with -mt, without loading the whole
    table, so a malformed metadata file stops the run before the analysis starts.

    Keyword arguments:
    metadata_file: csv table with taxonomy and other metadata, parsed with -mt parameter
    """
    header = pd.read_csv(metadata_file, nrows=0).columns.tolist()

    return check_metadata_columns(header)


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
