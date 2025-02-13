"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI saccharolipid (Lipids that contain a carbohydrate moiety.)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid that contains a carbohydrate moiety.
    Here we use two heuristics:
      1. Lipid part: the molecule must contain at least one long aliphatic chain 
         (here defined as a chain of 7 or more contiguous carbon atoms).
      2. Carbohydrate part: the molecule must contain at least one ring of 5 or 6 members
         that is composed solely of carbon and oxygen atoms (with at least one oxygen).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a saccharolipid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Detect lipid chain: a contiguous chain of 7 or more aliphatic carbons ---
    # We define a simple SMARTS for a chain of seven carbons.
    lipid_smarts = Chem.MolFromSmarts("CCCCCCC")
    lipid_matches = mol.GetSubstructMatches(lipid_smarts)
    if not lipid_matches:
        return False, "No long aliphatic (lipid) chain detected"
    
    # --- Detect carbohydrate moiety: search rings (5 or 6 members) with only C and O ---
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_found = False
    for ring in rings:
        if len(ring) in [5, 6]:
            # Check if each atom in the ring is either carbon (atomic number 6) or oxygen (atomic number 8)
            # and ensure that at least one oxygen is present.
            elems = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring]
            if all(el in (6, 8) for el in elems) and any(el == 8 for el in elems):
                sugar_found = True
                break

    if not sugar_found:
        return False, "No carbohydrate (sugar) moiety detected"

    # If both a lipid chain and a sugar ring are present, we classify it as a saccharolipid.
    return True, "Molecule contains both a long aliphatic (lipid) chain and a carbohydrate moiety"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided example SMILES strings for a saccharolipid 
    # (for instance, using the sugar-lipid structure from the given examples)
    smiles_example = "CCCCCCCCCCCCC(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(O)(=O)=O)[C@@H]1OC(=O)CCCCCCCCCCCCCCC"
    result, reason = is_saccharolipid(smiles_example)
    print("Is saccharolipid?", result)
    print("Reason:", reason)