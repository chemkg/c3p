"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3β-hydroxy-Δ(5)-steroid
A molecule is considered a 3β-hydroxy-Δ(5)-steroid if:
  1. It exhibits the fused steroid nucleus (heuristically: at least three 6-membered rings and one 5-membered ring).
  2. It contains a hydroxyl group on a ring carbon that is connected to at least two ring atoms (approximating the 3β-hydroxyl).
  3. It contains at least one double bond between two ring carbons (approximating the Δ(5) double bond).
Note: This is a heuristic approach.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy-Δ(5)-steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the steroid class
        str: Explanation for the classification decision
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a steroid-like fused ring system.
    # Heuristic: expect at least three six-membered rings and one five-membered ring.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    six_membered = sum(1 for ring in rings if len(ring) == 6)
    five_membered = sum(1 for ring in rings if len(ring) == 5)
    if six_membered < 3 or five_membered < 1:
        return False, ("Molecule does not appear to have a typical steroid fused ring system "
                       "(expected at least three 6-membered rings and one 5-membered ring, "
                       f"found {six_membered} six-membered and {five_membered} five-membered)")

    # Look for a hydroxyl group on a ring carbon.
    # The SMARTS "[C;R][OX2H]" finds an oxygen (with two connections) bonded to 
    # a carbon that is in a ring. Then, to ensure it is likely the 3β-OH, we check that
    # the carbon is connected to at least two other ring atoms.
    hydroxyl_smarts = "[C;R][OX2H]"
    hydroxyl_query = Chem.MolFromSmarts(hydroxyl_smarts)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_query)
    
    valid_hydroxyl = False
    for match in hydroxyl_matches:
        # match[0] is the index of the carbon atom with the OH attached.
        carbon_atom = mol.GetAtomWithIdx(match[0])
        ring_neighbor_count = 0
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.IsInRing():
                ring_neighbor_count += 1
        # In a steroid, the C3 bearing the OH is typically connected to multiple ring atoms.
        if ring_neighbor_count >= 2:
            valid_hydroxyl = True
            break
    if not valid_hydroxyl:
        return False, "No valid 3β-hydroxyl group found on a ring carbon of the steroid nucleus"

    # Look for a double bond within the ring system,
    # representing the Δ(5) double bond. We search for any double bond 
    # connecting two ring carbons.
    double_bond_smarts = "[C;R]=[C;R]"
    double_bond_query = Chem.MolFromSmarts(double_bond_smarts)
    db_matches = mol.GetSubstructMatches(double_bond_query)
    if not db_matches:
        return False, "No ring double bond (Δ(5)) found connecting two ring carbons"

    # If all criteria are met, then the molecule is classified as a 3β-hydroxy-Δ(5)-steroid.
    return True, "Molecule contains a steroid nucleus with a 3β-hydroxyl group and a Δ(5) double bond"


# Example usage:
if __name__ == "__main__":
    # Test with cholesterol SMILES (commonly a 3β-hydroxy-Δ(5)-steroid)
    chol_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC=C2C[C@H](C1)O)[H])" \
                  "(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])C"
    result, reason = is_3beta_hydroxy_Delta_5__steroid(chol_smiles)
    print("Cholesterol classification:", result)
    print("Reason:", reason)