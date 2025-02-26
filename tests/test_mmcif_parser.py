"""Tests the PDB parser."""

from collections import Counter
import pytest
import pathlib
import gzip

import ampal

TEST_FILE_FOLDER = pathlib.Path(__file__).parent / "testing_files"


def count_atom_lines(file_path: pathlib.Path, is_gzipped: bool = False) -> int:
    """Counts the number of lines starting with "ATOM " in an mmCIF file.

    Parameters
    ----------
    file_path : pathlib.Path
        Path to the mmCIF file.
    is_gzipped : bool, optional
        Whether the file is gzipped, by default False.

    Returns
    -------
    int
        The number of lines starting with "ATOM ".
    """
    if is_gzipped:
        with gzip.open(file_path, "rt") as file:
            count = 0
            for line in file:
                if "ATOM" in line or "HETATM" in line:
                    count += 1
    else:
        with open(file_path, "r") as file:
            count = 0
            for line in file:
                if "ATOM" in line or "HETATM" in line:
                    count += 1
    return count


def test_2j58():
    file_path = TEST_FILE_FOLDER / "2j58_1.cif"

    number_of_atoms = count_atom_lines(file_path, is_gzipped=False)
    ac = ampal.load_mmcif_file(file_path, is_gzipped=False)
    assert (
        len(list(ac[0].get_atoms())) == number_of_atoms  # type: ignore
    ), "Number of atoms in cif file and AMPAL object should match"
    return


def test_af_g5eb01():
    file_path = TEST_FILE_FOLDER / "AF-G5EB01-F1-model_v4.cif.gz"

    number_of_atoms = count_atom_lines(file_path, is_gzipped=True)
    ac = ampal.load_mmcif_file(file_path, is_gzipped=True)
    assert (
        len(list(ac[0].get_atoms())) == number_of_atoms  # type: ignore
    ), "Number of atoms in cif file and AMPAL object should match"
    return


def test_compare_2ht0():
    mmcif_file_path = TEST_FILE_FOLDER / "2ht0.cif"
    pdb_file_path = TEST_FILE_FOLDER / "2ht0.pdb"

    cif_assembly = ampal.load_mmcif_file(mmcif_file_path, is_gzipped=False)[0]
    pdb_assembly = ampal.load_pdb(pdb_file_path)

    assert len(list(cif_assembly.get_atoms(inc_alt_states=True))) == len(
        list(pdb_assembly.get_atoms(inc_alt_states=True))
    ), "Should have the same number of atoms."
    assert len(cif_assembly) == len(pdb_assembly), "Should have same number of chains"
    for chain in cif_assembly:
        assert len(cif_assembly[chain.id]) == len(
            pdb_assembly[chain.id]
        ), "Should have same number of monomers in chains"
