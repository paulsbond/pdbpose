from argparse import ArgumentParser, Namespace
from pathlib import Path
import gemmi
from . import pdbe


def parse_args() -> Namespace:
    parser = ArgumentParser()
    parser.add_argument("uniprot")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reference = None
    for pdb_id, asym_id in sorted(pdbe.chains(args.uniprot)):
        print(pdb_id, asym_id)
        structure = pdbe.structure(pdb_id)
        trimmed = trim(structure, asym_id)
        if reference is None:
            reference = trimmed
        else:
            superpose(fixed=reference, movable=trimmed)
        write_pdb(trimmed, pdb_id, asym_id)


def trim(structure: gemmi.Structure, asym_id: str) -> gemmi.Structure:
    trimmed = gemmi.Structure()
    model = trimmed.find_or_add_model("1")
    chain = model.add_chain("A")
    chain.append_residues(structure[0].get_subchain(asym_id))
    return trimmed


def superpose(fixed: gemmi.Structure, movable: gemmi.Structure) -> None:
    polymer1 = fixed[0].subchains()[0]
    polymer2 = movable[0].subchains()[0]
    polymer_type = gemmi.PolymerType.PeptideL
    selection = gemmi.SupSelect.CaP
    result = gemmi.calculate_superposition(polymer1, polymer2, polymer_type, selection)
    for chain in movable[0]:
        for residue in chain:
            for atom in residue:
                atom.pos = gemmi.Position(result.transform.apply(atom.pos))


def write_pdb(structure: gemmi.Structure, pdb_id: str, asym_id: str) -> None:
    path = Path("pdbpose", f"{pdb_id}_{asym_id}.pdb")
    path.parent.mkdir(exist_ok=True)
    structure.write_minimal_pdb(str(path))
