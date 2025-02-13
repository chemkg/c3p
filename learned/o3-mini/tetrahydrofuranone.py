"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: tetrahydrofuranone (any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    
    A tetrahydrofuranone is defined here as any molecule containing a five-membered ring
    (an oxolane ring) that has exactly one oxygen and four carbons (allowing some modifications)
    with at least one of the ring carbon atoms bearing a carbonyl group (C=O) as a substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a tetrahydrofuranone, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information (list of rings with atom indices)
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    
    # Iterate over all rings in the molecule
    for ring in rings:
        if len(ring) == 5:
            # Count number of oxygen and carbon atoms in the ring.
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            carbon_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # We expect an oxolane ring to have exactly one oxygen and four carbons.
            if oxygen_count == 1 and carbon_count == 4:
                # Check if any ring carbon bears an oxo (-C=O) substituent.
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    # We only need to check carbon atoms.
                    if atom.GetAtomicNum() == 6:
                        # Iterate over neighbors of the ring carbon.
                        for nbr in atom.GetNeighbors():
                            # Skip if the neighbor is within the ring (unless desired to check in-ring carbonyl)
                            if nbr.GetIdx() in ring:
                                continue
                            # Check if the neighbor is oxygen.
                            if nbr.GetAtomicNum() == 8:
                                # Get the bond between the ring carbon and the oxygen neighbor.
                                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                                if bond is not None:
                                    # A double bond (C=O) usually has a bond order of 2.0.
                                    # (BondTypeAsDouble will give a float value indicating bond multiplicity.)
                                    if bond.GetBondTypeAsDouble() >= 2.0:
                                        return True, "Found a tetrahydrofuran ring with an oxo substituent on a ring carbon"
                # Special note: if no exocyclic carbonyl is found attached to a ring carbon, we continue checking other rings.
    return False, "No tetrahydrofuran ring with an oxo substituent found"
    
# Example usage:
if __name__ == "__main__":
    # Testing with one of the provided examples: N-isovaleryl-L-homoserine lactone
    test_smiles = "CC(C)CC(=O)N[C@H]1CCOC1=O"
    result, reason = is_tetrahydrofuranone(test_smiles)
    print(f"SMILES: {test_smiles}")
    print("Classification:", result)
    print("Reason:", reason)