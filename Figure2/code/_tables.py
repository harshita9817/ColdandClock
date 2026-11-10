"""Shared table I/O for the preparation and enrichment scripts.

Generated tables are written atomically and deterministically: the gzip member
carries ``mtime=0`` so that rerunning an unchanged analysis reproduces the same
bytes, and the final rename only happens once the whole file is on disk.
"""

from __future__ import annotations

import gzip
import io
import os
import tempfile
from pathlib import Path

import pandas as pd


def read_table(path: Path, separator: str = ",") -> pd.DataFrame:
    """Read a (optionally gzipped) delimited table with literal string values.

    Gene and term identifiers are kept as written: no dtype inference, and no
    conversion of strings such as ``NA`` into missing values.
    """
    return pd.read_csv(path, sep=separator, dtype=str, keep_default_na=False)


def write_table(path: Path, frame: pd.DataFrame, separator: str = "\t") -> None:
    """Write ``frame`` to ``path`` as a deterministic, atomically replaced gzip table."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as raw_handle:
            with gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as gzip_handle:
                with io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as text_handle:
                    frame.to_csv(text_handle, sep=separator, index=False, lineterminator="\n")
        os.replace(temporary_name, path)
        os.chmod(path, 0o664)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)
