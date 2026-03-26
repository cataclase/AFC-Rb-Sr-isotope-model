import pandas as pd
import numpy as np
import re
import unicodedata


def strip_accents(s):
    s = unicodedata.normalize("NFKD", str(s))
    return "".join(ch for ch in s if not unicodedata.combining(ch))


def normalize_series(x):
    if pd.isna(x):
        return np.nan

    s = strip_accents(str(x)).lower().strip()
    key = re.sub(r"[^a-z]", "", s)

    mapping = {
        "alk": "Alk",
        "alkaline": "Alk",
        "phkcalcalk": "PHKCalcAlk",
        "phkcalcalkalk": "PHKCalcAlk",
        "ehkcalcalk": "EHKCalcAlk",
        "ehkcalcalkalk": "EHKCalcAlk",
        "calcalk": "CalcAlk",
        "chalk": "ChAlk",
        "shos": "Shos",
        "shoshonitic": "Shos",
        "archean": "Archean",
        "arqueano": "Archean",
        "paleoproterozoic": "Paleoproterozoic",
        "paleoproterozoico": "Paleoproterozoic",
    }

    return mapping.get(key, str(x).strip())


def clean_numeric(x, below_detection="value"):
    """
    below_detection:
        - 'value' -> <0.10 vira 0.10
        - 'half'  -> <0.10 vira 0.05
        - 'nan'   -> <0.10 vira NaN
    """
    if pd.isna(x):
        return np.nan

    s = str(x).strip()

    if s == "":
        return np.nan

    s = s.replace("−", "-").replace(",", ".")

    if s.lower() in {"bdl", "n.d.", "nd", "na", "nan", "none"}:
        return np.nan

    if s.startswith("<"):
        raw = re.sub(r"[^0-9eE\+\-\.]", "", s[1:])
        if raw == "":
            return np.nan
        val = float(raw)
        if below_detection == "half":
            return val / 2
        if below_detection == "nan":
            return np.nan
        return val

    s = re.sub(r"[^0-9eE\+\-\.]", "", s)
    if s == "":
        return np.nan

    return float(s)


def clean_numeric_series(series, below_detection="value"):
    return series.apply(lambda x: clean_numeric(x, below_detection=below_detection))