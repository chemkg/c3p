"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python3
"""
Classifies: 7-hydroxyisoflavones
Definition: A hydroxyisoflavone compound having a hydroxy (–OH) group at the 7–position (i.e. on the A ring of the isoflavone core).
The algorithm follows these steps:
 1. Verify that the molecule contains an isoflavone core.
    For this we use a SMARTS pattern for a 3–phenylchromen-4-one, which is a reasonably specific representation.
 2. Among the aromatic six-membered rings in the molecule, identify those that are fused (i.e. they share atoms with at least one other ring).
    Then check if at least one atom in such a ring carries an –OH group.
If both conditions are met, the molecule is classified as a 7-hydroxyisoflavone.
Note: This heuristic procedure has limitations and may not correctly classify all edge cases.
"""

from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    It first checks for the presence of an isoflavone core (3-phenylchromen-4-one)
    and then verifies that one of the aromatic rings that is fused within the core carries an -OH substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a 7-hydroxyisoflavone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS for the isoflavone core.
    # This pattern represents a 3–phenylchromen-4-one (isoflavone) core.
    # The SMARTS: "c1ccc2c(c1)oc(=O)c3ccccc23" will match a fused bicyclic system (chromen-4-one) with an attached phenyl ring.
    iso_core_smarts = "c1ccc2c(c1)oc(=O)c3ccccc23"
    iso_core = Chem.MolFromSmarts(iso_core_smarts)
    if iso_core is None:
        return False, "Error in SMARTS pattern for isoflavone core."
    
    if not mol.HasSubstructMatch(iso_core):
        return False, "Molecule does not contain the expected isoflavone core."
    
    # Pre-calculate ring membership counts for all atoms.
    ring_info = mol.GetRingInfo().AtomRings()
    atom_ring_counts = {atom_idx:0 for atom_idx in range(mol.GetNumAtoms())}
    for ring in ring_info:
        for atom_idx in ring:
            atom_ring_counts[atom_idx] += 1

    # Now search through the rings of the molecule.
    # We focus on aromatic, six-membered rings that are fused (i.e. some atoms are in more than one ring).
    for ring in ring_info:
        if len(ring) != 6:
            continue  # Only consider six–membered rings.
        # Check that all atoms in the ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Check if the ring is fused with another ring: at least one atom in the ring belongs to more than one ring.
        if not any(atom_ring_counts[idx] > 1 for idx in ring):
            continue  # Likely an isolated phenyl (B ring), not part of the core.
        
        # For each atom in the ring, check for a hydroxyl (-OH) substituent.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # Expect carbon.
                continue
            # Look at neighbors: an -OH group is an oxygen with at least one attached hydrogen.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    return True, "Molecule contains an isoflavone core with an -OH on the fused aromatic (A) ring."
    
    # If no fused aromatic ring carries an -OH:
    return False, "Molecule contains an isoflavone core but no hydroxyl group was found on a fused (A ring) benzene."

# Uncomment the lines below for a quick test.
# if __name__ == "__main__":
#     # Example: 7-hydroxyisoflavone
#     test_smiles = "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O"
#     result, reason = is_7_hydroxyisoflavones(test_smiles)
#     print(result, reason)