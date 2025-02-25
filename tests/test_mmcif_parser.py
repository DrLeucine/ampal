"""Tests the PDB parser."""

from collections import Counter
import pytest
import pathlib

import ampal

TEST_FILE_FOLDER = pathlib.Path(__file__).parent / "testing_files"


def test_2j52():
    file_path = TEST_FILE_FOLDER / "2j58_1.cif"
    ampal.load_mmcif_file(file_path, is_gzipped=False)
    return


def test_af_g5eb01():
    file_path = TEST_FILE_FOLDER / "AF-G5EB01-F1-model_v4.cif.gz"
    ampal.load_mmcif_file(file_path, is_gzipped=True)
    return
