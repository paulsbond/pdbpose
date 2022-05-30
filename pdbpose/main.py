from argparse import ArgumentParser, Namespace
from pathlib import Path
import gemmi
from . import pdbe


def parse_args() -> Namespace:
    parser = ArgumentParser()
    parser.add_argument("uniprot")
    parser.add_argument("--radius", type=float, default=4.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reference = None
    for pdb_id, asym_id in sorted(pdbe.subchains(args.uniprot)):
        print(pdb_id, asym_id)
        structure = pdbe.structure(pdb_id)
        trimmed = trim(structure, asym_id, args.radius)
        if reference is None:
            reference = trimmed
        else:
            superpose(fixed=reference, movable=trimmed)
        write_pdb(trimmed, pdb_id, asym_id)


def trim(structure: gemmi.Structure, asym_id: str, radius: float) -> gemmi.Structure:
    trimmed = gemmi.Structure()
    model = trimmed.find_or_add_model("1")
    model.add_chain("A")
    model.add_chain("B")
    search = gemmi.NeighborSearch(structure, max_radius=radius)
    search.populate(include_h=False)
    transforms = {}
    for residue in structure[0].get_subchain(asym_id):
        model["A"].add_residue(residue)
        for atom in residue:
            for mark in search.find_neighbors(atom, max_dist=radius):
                cra = mark.to_cra(structure[0])
                transform = search.get_image_transformation(mark.image_idx)
                image = structure.cell.find_nearest_pbc_image(
                    atom.pos, cra.atom.pos, mark.image_idx
                )
                vec = transform.vec + gemmi.Vec3(*image.pbc_shift)
                transform = gemmi.Transform(transform.mat, vec)
                if cra.residue.subchain == asym_id and transform.is_identity():
                    continue
                added = transforms.setdefault((mark.chain_idx, mark.residue_idx), [])
                if any(t.approx(transform, epsilon=1e-5) for t in added):
                    continue
                added.append(transform)
                cloned_residue = cra.residue.clone()
                num = model["B"][-1].seqid.num + 1 if len(model["B"]) > 0 else 1
                cloned_residue.seqid = gemmi.SeqId(num, " ")
                for cloned_atom in cloned_residue:
                    fractional = structure.cell.fractionalize(cloned_atom.pos)
                    transformed = gemmi.Fractional(*transform.apply(fractional))
                    cloned_atom.pos = structure.cell.orthogonalize(transformed)
                model["B"].add_residue(cloned_residue)
    trimmed.remove_empty_chains()
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
