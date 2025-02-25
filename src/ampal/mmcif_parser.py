"""Contains code for parsing mmCIF files."""

from collections import OrderedDict
import itertools
from typing import Tuple, List, Dict
import gzip
import pathlib

from .base_ampal import Atom, Monomer
from .assembly import AmpalContainer, Assembly
from .protein import Polypeptide, Residue
from .nucleic_acid import Polynucleotide, Nucleotide
from .ligands import Ligand, LigandGroup
from .amino_acids import standard_amino_acids
from .data import PDB_ATOM_COLUMNS


def load_mmcif_file(
    file_path: pathlib.Path, is_gzipped: bool = False
) -> AmpalContainer:
    if is_gzipped:
        file = gzip.open(file_path)
    else:
        file = open(file_path, "r")

    parsing = False
    column_labels = []
    atom_lines = []
    for in_line in file.readlines():
        line = str(in_line)
        if not parsing:
            if line.startswith("_atom_site."):
                parsing = True
                column_labels.append(line.split(".")[1].strip())
        else:
            if line.startswith("_atom_site."):
                column_labels.append(line.split(".")[1].strip())
            elif line.startswith("#"):
                parsing = False
            else:
                atom_lines.append(line)

    # at_type = line[0:6].strip()  # 0
    # at_ser = int(line[6:11].strip())  # 1
    # at_name = line[12:16].strip()  # 2
    # alt_loc = line[16].strip()  # 3
    # res_name = line[17:20].strip()  # 4
    # chain_id = line[21].strip()  # 5
    # res_seq = int(line[22:26].strip())  # 6
    # i_code = line[26].strip()  # 7
    # x = float(line[30:38].strip())  # 8
    # y = float(line[38:46].strip())  # 9
    # z = float(line[46:54].strip())  # 10
    # occupancy = float(line[54:60].strip())  # 11
    # _temp_factor = line[60:66].strip()
    # if _temp_factor:
    #     temp_factor = float(_temp_factor)
    # else:
    #     temp_factor = 0.0
    # element = line[76:78].strip()  # 13
    # charge = line[78:80].strip()  # 14
    states = {}

    record_type_index = column_labels.index("group_PDB")
    at_ser_index = column_labels.index("id")
    at_name_index = column_labels.index("label_atom_id")
    alt_loc_index = column_labels.index("label_alt_id")
    res_name_index = column_labels.index("label_comp_id")
    chain_id_index = column_labels.index("label_asym_id")
    res_seq_id_index = column_labels.index("auth_seq_id")
    i_code_index = column_labels.index("pdbx_PDB_ins_code")
    x_index = column_labels.index("Cartn_x")
    y_index = column_labels.index("Cartn_y")
    z_index = column_labels.index("Cartn_z")
    occupancy_index = column_labels.index("occupancy")
    bfactor_index = column_labels.index("B_iso_or_equiv")
    element_index = column_labels.index("type_symbol")
    charge_index = column_labels.index("pdbx_formal_charge")
    model_number_index = column_labels.index("pdbx_PDB_model_num")

    for atom_line in atom_lines:
        atom_columns = atom_line.strip().split(" ")
        record_type = atom_columns[record_type_index]
        model_number = atom_columns[model_number_index]
        if model_number not in states:
            states[model_number] = {}
        state = states[model_number]

        chain_id = atom_columns[chain_id_index]
        if chain_id not in state:
            state[chain_id] = ({}, set())
        (chain, chain_composition) = state[chain_id]

        res_seq_id = atom_columns[res_seq_id_index]
        i_code = atom_columns[i_code_index]
        full_res_id = (res_seq_id, i_code)
        if full_res_id not in chain:
            chain[full_res_id] = {}
        residue = chain[full_res_id]

        res_name = atom_columns[res_name_index]
        res_label = atom_columns[at_name_index]
        atom_coordinate = (
            float(atom_columns[x_index]),
            float(atom_columns[y_index]),
            float(atom_columns[z_index]),
        )
        atom = Atom(
            atom_coordinate,
            element=atom_columns[element_index],
            atom_id=atom_columns[at_ser_index],
            res_label=res_label,
            occupancy=atom_columns[occupancy_index],
            bfactor=float(atom_columns[bfactor_index]),
            charge=atom_columns[charge_index],
            state=atom_columns[i_code_index],
            parent=None,
        )
        assert (
            res_name,
            res_label,
        ) not in residue, "Atom label is not unique, going to overwrite!"

        residue[(res_name, res_label)] = atom
        if record_type == "ATOM":
            if res_name in standard_amino_acids.values():
                chain_composition.add("P")
            else:
                chain_composition.add("N")
        else:
            chain_composition.add("H")

    ampal_container = AmpalContainer(id=str(file_path.stem))
    for state_id, state in states.items():
        assembly = Assembly(assembly_id=state_id)
        for chain_id, (chain, chain_composition) in state.items():
            if "P" in chain_composition:
                polymer = Polypeptide(polymer_id=chain_id, parent=assembly)
                polymer.ligands = LigandGroup()
            elif "N" in chain_composition:
                polymer = Polynucleotide(polymer_id=chain_id, parent=assembly)
                polymer.ligands = LigandGroup()
            else:
                polymer = LigandGroup(polymer_id=chain_id, parent=assembly)
            for (res_seq_id, i_code), residue in chain.items():
                residue_labels = list({x[0] for x in residue.keys()})
                assert (
                    len(residue_labels) == 1
                ), "A single residue has has multiple residue labels"
                mol_code = residue_labels[0]
                if isinstance(polymer, Polypeptide):
                    if mol_code in standard_amino_acids.values():
                        monomer = Residue(
                            monomer_id=res_seq_id,
                            mol_code=mol_code,
                            insertion_code=i_code,
                            is_hetero=False,
                            parent=polymer,
                        )
                        polymer.append(monomer)
                    else:
                        monomer = Ligand(
                            monomer_id=res_seq_id,
                            mol_code=mol_code,
                            insertion_code=i_code,
                            is_hetero=True,
                            parent=polymer,
                        )
                        polymer.ligands.append(monomer)
                elif isinstance(polymer, Polynucleotide):
                    monomer = Nucleotide(
                        monomer_id=res_seq_id,
                        mol_code=mol_code,
                        insertion_code=i_code,
                        is_hetero=False,
                        parent=polymer,
                    )
                    if mol_code in ["U", "G", "C", "A", "DT", "DG", "DC", "DA"]:
                        polymer.append(monomer)
                    else:
                        polymer.ligands.append(monomer)

                else:
                    monomer = Ligand(
                        monomer_id=res_seq_id,
                        mol_code=mol_code,
                        insertion_code=i_code,
                        is_hetero=True,
                        parent=polymer,
                    )
                    polymer.append(monomer)
                for (_, atom_label), atom in residue.items():
                    atom.parent = monomer
                    monomer[atom_label] = atom
            assembly.append(polymer)
        ampal_container.append(assembly)

    print(ampal_container)
    print(ampal_container[0])
    print(ampal_container[0]["E"])
    print(len(list(ampal_container[0].get_atoms(ligands=True, inc_alt_states=True))))

    file.close()
    return ampal_container
