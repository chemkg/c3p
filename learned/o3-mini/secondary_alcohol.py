"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
Definition: A secondary alcohol is a compound in which a hydroxy group (-OH) 
is attached to a saturated carbon atom (sp3) that has two carbon substituents and one hydrogen.

This function uses RDKit to parse a SMILES string and check for the secondary alcohol pattern.
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule contains a secondary alcohol group.
    
    A secondary alcohol group is defined as an sp3 carbon bonded to:
      - exactly one hydroxyl group (-OH) [oxygen bonded to one hydrogen],
      - exactly two carbon atoms (R groups),
      - exactly one hydrogen.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if at least one secondary alcohol group is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into a molecule using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that hydrogen counts are accurate.
    mol = Chem.AddHs(mol)

    # Iterate over all atoms to search for a candidate carbon bearing an -OH group.
    for atom in mol.GetAtoms():
        # Check if the atom is a carbon and is tetrahedral (sp3).
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            # We will count the different types of neighbors for this carbon.
            oxygen_count = 0
            carbon_count = 0
            hydrogen_count = 0

            # Loop over neighbors of the carbon atom.
            for neighbor in atom.GetNeighbors():
                atomic_num = neighbor.GetAtomicNum()
                if atomic_num == 8:
                    # Check if this oxygen is part of a hydroxyl group:
                    # In an -OH group, the oxygen should have exactly one hydrogen neighbor.
                    h_on_oxygen = sum(1 for nb in neighbor.GetNeighbors() if nb.GetAtomicNum() == 1)
                    if h_on_oxygen == 1:
                        oxygen_count += 1
                elif atomic_num == 6:
                    carbon_count += 1
                elif atomic_num == 1:
                    hydrogen_count += 1

            # For a secondary alcohol at this carbon, we expect:
            # one -OH neighbor, two carbon neighbors, and one hydrogen.
            if oxygen_count == 1 and carbon_count == 2 and hydrogen_count == 1:
                return True, f"Molecule contains a secondary alcohol group at carbon atom index {atom.GetIdx()}."

    return False, "No secondary alcohol group found."