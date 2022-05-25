from argparse import ArgumentParser, Namespace
from pathlib import Path
from gemmi import Structure
from . import pdbe


def parse_args() -> Namespace:
    parser = ArgumentParser()
    parser.add_argument("uniprot")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for pdb_id, asym_id in sorted(pdbe.chains(args.uniprot)):
        print(pdb_id, asym_id)
        structure = pdbe.structure(pdb_id)
        trimmed = trim(structure, asym_id)
        write_pdb(trimmed, pdb_id, asym_id)


def trim(structure: Structure, asym_id: str) -> Structure:
    trimmed = Structure()
    model = trimmed.find_or_add_model("1")
    chain = model.add_chain("A")
    chain.append_residues(structure[0].get_subchain(asym_id))
    return trimmed


def write_pdb(structure: Structure, pdb_id: str, asym_id: str) -> None:
    path = Path("pdbpose", f"{pdb_id}_{asym_id}.pdb")
    path.parent.mkdir(exist_ok=True)
    structure.write_minimal_pdb(str(path))
