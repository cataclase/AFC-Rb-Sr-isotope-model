from pathlib import Path
import pandas as pd


DEFAULT_EXAMPLE_PATH = Path("examples") / "Example_DRN.xlsx"
DEFAULT_SHEET = "DRN_magmatism"


def load_example_dataset(path=DEFAULT_EXAMPLE_PATH, sheet_name=DEFAULT_SHEET):
    return pd.read_excel(path, sheet_name=sheet_name)


def load_geochemistry(path=DEFAULT_EXAMPLE_PATH, sheet_name=DEFAULT_SHEET):
    return pd.read_excel(path, sheet_name=sheet_name)


def load_rb_sr(path=DEFAULT_EXAMPLE_PATH, sheet_name=DEFAULT_SHEET):
    return pd.read_excel(path, sheet_name=sheet_name)


def load_sm_nd(path=DEFAULT_EXAMPLE_PATH, sheet_name=DEFAULT_SHEET):
    return pd.read_excel(path, sheet_name=sheet_name)