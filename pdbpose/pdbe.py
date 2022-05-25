from pathlib import Path
from typing import Iterable, Tuple
import gemmi
import requests


SERVER = "https://www.ebi.ac.uk/pdbe"


def chains(uniprot: str) -> Iterable[Tuple[str, str]]:
    url = f"{SERVER}/search/pdb/select?"
    query = f"uniprot_accession:{uniprot} AND status:REL"
    filter_list = "pdb_id,struct_asym_id"
    data = {"q": query, "fl": filter_list, "rows": 1000000, "wt": "json"}
    response = requests.post(url, data=data, timeout=30)
    response.raise_for_status()
    for doc in response.json().get("response", {}).get("docs", []):
        for asym_id in doc["struct_asym_id"]:
            yield doc["pdb_id"], asym_id


def structure(pdb_id: str) -> gemmi.Structure:
    filename = f"{pdb_id}.cif"
    path = Path("downloads", filename)
    if not path.exists():
        path.parent.mkdir(exist_ok=True)
        url = f"{SERVER}/entry-files/download/{filename}"
        response = requests.get(url, timeout=30, stream=True)
        response.raise_for_status()
        with open(path, "wb") as stream:
            for chunk in response.iter_content():
                stream.write(chunk)
    return gemmi.read_structure(str(path))
