import gzip
import pathlib
from collections import OrderedDict
from typing import Dict, List, Tuple, TextIO, Set, TypedDict, Union
import warnings

from ampal.base_ampal import Atom
from ampal.assembly import AmpalContainer, Assembly
from ampal.protein import Polypeptide, Residue
from ampal.nucleic_acid import Polynucleotide, Nucleotide
from ampal.ligands import Ligand, LigandGroup
from ampal.amino_acids import standard_amino_acids

# Type aliases for improved readability and maintainability
ChainDict = Dict[Tuple[str, str], Dict[Tuple[str, str, str], Atom]]
ChainComposition = Set[str]  # Contains 'P', 'N', or 'H'
StateDict = Dict[str, Tuple[ChainDict, ChainComposition]]
StatesDict = Dict[str, StateDict]


# TypedDict for atom data structure, enhancing type safety
class AtomData(TypedDict):
    """
    Represents the nested structure of atom data extracted from an mmCIF file.

    Attributes
    ----------
    model_number : str
        Model number in which chain is present
    chain_id : str
        A dictionary representing a chain, containing residues and chain type.
        The keys are chain identifiers (typically single letters).  The values
        are tuples. The first element of the tuple is a dictionary of Residues
        keyed by residue ID. The second element is a set representing the
        chain composition (containing 'P' for protein, 'N' for nucleic acid,
        'H' for hetero).
    res_seq_id : str
        Unique identifier for residues that contains the residue number and the insertion code
    residue : str
        A dictionary containing all atoms for a single Residue.
        The keys are tuples: (residue name, atom label). The values are Atom
        objects.
    """

    model_number: str
    chain_id: str
    res_seq_id: Tuple[str, str]
    residue: Dict[Tuple[str, str], Atom]


def _open_mmcif_file(file_path: pathlib.Path, is_gzipped: bool = False) -> TextIO:
    """Opens an mmCIF file, handling gzipped files if necessary.

    Parameters
    ----------
    file_path : pathlib.Path
        Path to the mmCIF file.
    is_gzipped : bool, optional
        Whether the file is gzipped, by default False.

    Returns
    -------
    TextIO
        File object.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    try:
        if is_gzipped:
            return gzip.open(file_path, "rt")  # Open in text mode
        else:
            return open(file_path, "r")
    except OSError as e:
        raise OSError(f"Error opening file {file_path}: {e}") from e
    except Exception as e:
        raise Exception(
            f"An unexpected error occurred while opening {file_path}: {e}"
        ) from e


def _parse_atom_site_records(file: TextIO) -> Tuple[List[str], List[str]]:
    """Parses the _atom_site records in an mmCIF file.

    Parameters
    ----------
    file : TextIO
        File object of the mmCIF file.

    Returns
    -------
    Tuple[List[str], List[str]]
        A tuple containing:
            - List of column labels.
            - List of atom lines.
    """
    parsing = False
    column_labels: List[str] = []
    atom_lines: List[str] = []
    for line in file:
        line = str(line)
        if not parsing:
            if line.startswith("_atom_site."):
                parsing = True
                column_labels.append(line.split(".")[1].strip())
        else:
            if line.startswith("_atom_site."):
                column_labels.append(line.split(".")[1].strip())
            elif line.startswith("#"):
                parsing = False  # Stop parsing at the end of the block
            elif (
                line.strip()
            ):  # Only add non-empty lines, avoids issues with whitespace
                atom_lines.append(line.strip())
    return column_labels, atom_lines


def _extract_atom_data(atom_lines: List[str], column_labels: List[str]) -> StatesDict:
    """Extracts atom data from atom lines and organizes it by model, chain, and residue.

    Parameters
    ----------
    atom_lines : List[str]
        List of atom lines.
    column_labels : List[str]
        List of column labels.

    Returns
    -------
    StatesDict
        A dictionary representing the structure:
        {model_number: {chain_id: ({ (res_seq_id, i_code): { (res_name, res_label): Atom } }, chain_composition)}}.  The chain composition is a set containing P, N or H.
    """
    # Use a try-except block to handle potential errors during index lookup
    try:
        indices = {
            label: column_labels.index(label)
            for label in [
                "group_PDB",
                "id",
                "label_atom_id",
                "label_alt_id",
                "label_comp_id",
                "auth_asym_id",
                "auth_seq_id",
                "pdbx_PDB_ins_code",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "occupancy",
                "B_iso_or_equiv",
                "type_symbol",
                "pdbx_formal_charge",
                "pdbx_PDB_model_num",
            ]
        }
    except ValueError as e:
        raise ValueError(f"Required column label not found in mmCIF file: {e}") from e

    states: StatesDict = {}
    for atom_line in atom_lines:
        # Use .split() without arguments to handle variable whitespace
        atom_columns = atom_line.split()
        # Remove quotes if present, handle missing columns robustly
        atom_columns = [col.replace('"', "") for col in atom_columns]

        # Helper function to get values with default handling
        def get_value(key, default=""):
            try:
                return atom_columns[indices[key]]
            except (IndexError, KeyError):
                return default

        record_type = get_value("group_PDB")
        model_number = get_value("pdbx_PDB_model_num")

        if model_number not in states:
            states[model_number] = {}
        state = states[model_number]

        chain_id = get_value("auth_asym_id")
        if chain_id not in state:
            state[chain_id] = ({}, set())
        chain, chain_composition = state[chain_id]

        res_seq_id = get_value("auth_seq_id")
        i_code = get_value("pdbx_PDB_ins_code")
        # Handle cases where i_code might be "?" or "." and should be treated as empty
        i_code = "" if i_code in ["?", "."] else i_code
        full_res_id = (res_seq_id, i_code)
        if full_res_id not in chain:
            chain[full_res_id] = {}
        residue = chain[full_res_id]

        res_name = get_value("label_comp_id")
        res_label = get_value("label_atom_id")
        alt_loc = get_value("label_alt_id", "")  # Get alt_loc, default to empty string
        alt_loc = "A" if alt_loc in ["?", "."] else alt_loc

        # Handle coordinate extraction with error checking
        try:
            atom_coordinate = (
                float(get_value("Cartn_x")),
                float(get_value("Cartn_y")),
                float(get_value("Cartn_z")),
            )
        except (ValueError, KeyError):
            warnings.warn(
                f"Invalid or missing coordinates for atom {res_label} in residue {res_seq_id}{i_code}, chain {chain_id}.  Skipping atom."
            )
            continue  # Skip this atom

        # Handle potential errors with charge and occupancy (e.g., if they are '?')
        charge = get_value("pdbx_formal_charge")
        charge = "" if charge in ["?", "."] else charge
        try:
            occupancy = float(get_value("occupancy", "1.0"))  # Default occupancy to 1.0
        except (ValueError, KeyError):
            occupancy = 1.0  # Default occupancy

        try:
            bfactor = float(get_value("B_iso_or_equiv"))
        except (ValueError, KeyError):
            warnings.warn(
                f"Invalid or missing B-factor for atom {res_label} in residue {res_seq_id}{i_code}, chain {chain_id}.  Setting to 0.0."
            )
            bfactor = 0.0

        try:
            atom = Atom(
                atom_coordinate,
                element=get_value("type_symbol"),
                atom_id=get_value("id"),
                res_label=res_label,
                occupancy=occupancy,
                bfactor=bfactor,
                charge=charge,
                state=alt_loc,  # Use alt_loc here
                parent=None,
            )
        except ValueError as e:
            warnings.warn(
                f"Error creating Atom object for {res_label} in residue {res_seq_id}{i_code}, chain {chain_id}: {e}.  Skipping atom."
            )
            continue

        if (res_name, res_label, alt_loc) in residue:
            warnings.warn(
                f"Atom label {(res_name, res_label, alt_loc)} is not unique, overwriting! {atom_line}"
            )

        residue[(res_name, res_label, alt_loc)] = atom
        if record_type == "ATOM":
            if res_name in standard_amino_acids.values():
                chain_composition.add("P")
            # More robust nucleic acid check
            elif res_name in [
                "A",
                "C",
                "G",
                "T",
                "U",
                "DA",
                "DC",
                "DG",
                "DT",
                "DU",
                "DI",
            ]:
                chain_composition.add("N")
            else:
                chain_composition.add("H")  # Assuming HETATM if not standard AA
        elif record_type == "HETATM":
            chain_composition.add("H")
        else:
            warnings.warn(
                f"Unknown record type: {record_type}.  Treating as hetero atom."
            )
            chain_composition.add("H")

    return states


def _create_ampal_structure(
    states: StatesDict, file_path: pathlib.Path
) -> AmpalContainer:
    """Creates an AmpalContainer structure from the extracted atom data.

    Parameters
    ----------
    states : StatesDict
        Dictionary containing atom data, as returned by _extract_atom_data.
    file_path : pathlib.Path
        Path of the file, to extract the name of the file to save as the name of the assembly.

    Returns
    -------
    AmpalContainer
        An AmpalContainer object representing the structure.
    """
    ampal_container = AmpalContainer(id=str(file_path.stem))
    for state_id, state in states.items():
        assembly = Assembly(assembly_id=f"{str(file_path.stem)}_{state_id}")
        for chain_id, (chain, chain_composition) in state.items():
            if "P" in chain_composition:
                polymer = Polypeptide(polymer_id=chain_id, parent=assembly)
                polymer.ligands = LigandGroup()  # Initialize LigandGroup
            elif "N" in chain_composition:
                polymer = Polynucleotide(polymer_id=chain_id, parent=assembly)
                polymer.ligands = LigandGroup()
            # Handle case with only hetero atoms: treat as a single LigandGroup
            else:  # "H" in chain_composition or empty
                polymer = LigandGroup(polymer_id=chain_id, parent=assembly)  # type: ignore

            for (res_seq_id, i_code), residue in chain.items():
                residue_labels = list({x[0] for x in residue.keys()})
                if len(residue_labels) != 1:
                    warnings.warn(
                        f"Residue {res_seq_id}{i_code} has multiple residue labels: {residue_labels}. Using the first one ({residue_labels[0]})."
                    )

                mol_code = residue_labels[0]

                if isinstance(polymer, (Polypeptide, Polynucleotide)):
                    if (
                        isinstance(polymer, Polypeptide)
                        and mol_code in standard_amino_acids.values()
                    ):
                        monomer = Residue(
                            monomer_id=res_seq_id,
                            mol_code=mol_code,
                            insertion_code=i_code,
                            is_hetero=False,
                            parent=polymer,
                        )
                        polymer.append(monomer)
                    elif isinstance(polymer, Polynucleotide) and mol_code in [
                        "U",
                        "G",
                        "C",
                        "A",
                        "DT",
                        "DG",
                        "DC",
                        "DA",
                        "DU",
                        "DI",
                    ]:
                        monomer = Nucleotide(
                            monomer_id=res_seq_id,
                            mol_code=mol_code,
                            insertion_code=i_code,
                            is_hetero=False,
                            parent=polymer,
                        )
                        polymer.append(monomer)

                    else:  # It is a ligand associated to a polymer
                        monomer = Ligand(
                            monomer_id=res_seq_id,
                            mol_code=mol_code,
                            insertion_code=i_code,
                            is_hetero=True,
                            parent=polymer,
                        )
                        if isinstance(polymer, (Polypeptide, Polynucleotide)):
                            polymer.ligands.append(monomer)  # type: ignore
                        else:
                            # Should not happen based on current logic, but included for completeness
                            polymer.append(monomer)

                else:  # polymer is a LigandGroup
                    monomer = Ligand(
                        monomer_id=res_seq_id,
                        mol_code=mol_code,
                        insertion_code=i_code,
                        is_hetero=True,
                        parent=polymer,
                    )
                    polymer.append(monomer)

                monomer.states = gen_states(list(residue.items()))
            assembly._molecules.append(polymer)
        ampal_container.append(assembly)

    return ampal_container


def gen_states(atoms: List[Tuple[Tuple[str, str, str], Atom]]) -> OrderedDict:
    """Generates the `states` dictionary for a `Monomer`.

    atoms : [Atom]
        A list of atom data parsed from the input PDB.
    """
    states = OrderedDict()
    for (_, atom_label, a_state_label), atom in atoms:
        state_label = "A" if not a_state_label else a_state_label
        if state_label not in states:
            states[state_label] = OrderedDict()
        states[state_label][atom_label] = atom

    # This code is to check if there are alternate states and populate any
    # both states with the full complement of atoms
    states_len = [(k, len(x)) for k, x in states.items()]
    if (len(states) > 1) and (len(set([x[1] for x in states_len])) > 1):
        for t_state, t_state_d in states.items():
            new_s_dict = OrderedDict()
            for k, v in states[sorted(states_len, key=lambda x: x[0])[0][0]].items():
                if k not in t_state_d:
                    c_atom = Atom(
                        v._vector,
                        v.element,
                        atom_id=v.id,
                        res_label=v.res_label,
                        occupancy=v.tags["occupancy"],
                        bfactor=v.tags["bfactor"],
                        charge=v.tags["charge"],
                        state=t_state[0],
                        parent=v.parent,
                    )
                    new_s_dict[k] = c_atom
                else:
                    new_s_dict[k] = t_state_d[k]
            states[t_state] = new_s_dict
    return states


def load_mmcif_file(
    file_path: Union[str, pathlib.Path], is_gzipped: bool = False
) -> AmpalContainer:
    """Loads and parses an mmCIF file, returning an AmpalContainer.

    Parameters
    ----------
    file_path : Union[str, pathlib.Path]
        Path to the mmCIF file, either as a string or a pathlib.Path object.
    is_gzipped : bool, optional
        Whether the file is gzipped, by default False.

    Returns
    -------
    AmpalContainer
        An AmpalContainer object representing the structure.

    Raises
    ------
    ValueError
        If the file format is invalid.
    TypeError
        If file_path is not a string or pathlib.Path.
    """
    if isinstance(file_path, str):
        file_path_obj = pathlib.Path(file_path)
    elif isinstance(file_path, pathlib.Path):
        file_path_obj = file_path
    else:
        raise TypeError("file_path must be a string or pathlib.Path object")

    try:
        with _open_mmcif_file(file_path_obj, is_gzipped) as file:
            column_labels, atom_lines = _parse_atom_site_records(file)
            if not column_labels or not atom_lines:
                raise ValueError("Invalid mmCIF file format: no atom data found.")
            states = _extract_atom_data(atom_lines, column_labels)
        ampal_container = _create_ampal_structure(states, file_path_obj)
    except Exception as e:
        raise ValueError(f"Error processing mmCIF file {file_path_obj}: {e}") from e

    return ampal_container


__author__ = "Christopher W. Wood"
