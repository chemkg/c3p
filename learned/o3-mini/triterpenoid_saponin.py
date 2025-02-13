"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: Triterpenoid saponin (a terpene glycoside in which the terpene moiety is a triterpenoid)

A triterpenoid saponin is defined as a molecule that contains one or more sugar (glycoside) units
attached via an O–glycosidic bond to a triterpenoid aglycone. Triterpenoid aglycones typically
have about 30 carbons and a fused ring system (usually 4 to 5 rings).

This implementation uses heuristics: 
  – It looks for one or more sugar rings, which are defined as 5- or 6-membered rings containing exactly one oxygen atom in the ring and multiple hydroxyl (–OH) substituents on the carbons.
  – It then “removes” the atoms from the sugar(s) and counts the number of carbon atoms in the remainder.
  – It also counts the ring systems exclusively composed of aglycone atoms.
If the aglycone appears to have roughly 30 carbons (we allow a little margin) and at least 4 fused rings,
we assume the molecule is a triterpenoid saponin.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    
    A triterpenoid saponin is defined as a terpene glycoside in which the terpene moiety is a triterpenoid.
    This heuristic checks for:
      - The presence of at least one sugar (glycoside) moiety.
      - An aglycone (non-sugar part) with about 30 carbon atoms.
      - A fused ring system (at least 4 rings) in the aglycone.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    # First, try to detect sugar rings.
    # We use a heuristic: a sugar ring is typically 5- or 6-membered, containing exactly one oxygen atom in the ring,
    # and its carbon atoms often bear one or more hydroxyl groups.
    sugar_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) in [5, 6]:
            # Count oxygen atoms present in the ring.
            num_oxy = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if num_oxy != 1:
                continue
            # Count hydroxyl substituents attached to carbons in the ring.
            hydroxyl_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # carbon atom
                    # Check neighbors: a hydroxyl will be an oxygen with degree 1.
                    for nbr in atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                            hydroxyl_count += 1
            if hydroxyl_count >= 2:
                sugar_rings.append(set(ring))
    
    if not sugar_rings:
        return False, "No sugar moiety (glycoside unit) detected"
    
    # Combine all atoms that belong to any detected sugar ring.
    sugar_atoms = set()
    for s in sugar_rings:
        sugar_atoms.update(s)
    
    # The aglycone is taken as the set of atoms not assigned to any sugar ring.
    total_atoms = set(range(mol.GetNumAtoms()))
    aglycone_atoms = total_atoms - sugar_atoms
    
    # Count carbon atoms in the aglycone.
    aglycone_carbons = sum(1 for idx in aglycone_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # Triterpenoids are built from 30 carbons (with a margin of ±3).
    if not (27 <= aglycone_carbons <= 33):
        return False, f"Aglycone carbon count ({aglycone_carbons}) not within expected range for a triterpenoid (~30 carbons)"
    
    # Count ring systems that are entirely within the aglycone.
    aglycone_ring_count = 0
    for ring in ring_info.AtomRings():
        if set(ring).issubset(aglycone_atoms):
            aglycone_ring_count += 1
    
    # Many triterpenoids are pentacyclic (5 rings) but allowing a bit of flexibility we require at least 4.
    if aglycone_ring_count < 4:
        return False, f"Fused ring system in aglycone insufficient (found {aglycone_ring_count} rings, expected at least 4)"
    
    return True, "Molecule contains a triterpenoid aglycone (~30 carbons with fused rings) and a sugar moiety"

# Example usage (you can remove or modify the following lines):
if __name__ == "__main__":
    # Example SMILES for ginsenoside Re (a triterpenoid saponin) provided in the prompt.
    test_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@](C)(CCC=C(C)C)O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_triterpenoid_saponin(test_smiles)
    print(result, ":", reason)