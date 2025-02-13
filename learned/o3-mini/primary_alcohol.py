"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: Primary Alcohol
Definition:
  A primary alcohol is defined as a compound in which a hydroxy (-OH) group is attached
  to a saturated carbon atom (sp3) that has either (a) three hydrogen atoms attached (CH3–OH)
  or (b) two hydrogen atoms and exactly one carbon neighbor (RCH2–OH).

This function attempts to detect at least one primary alcohol moiety in the molecule.
It does so by (1) converting the molecule to include all explicit hydrogens,
and (2) iterating over oxygen atoms that appear to be part of a simple –OH group.
Upon finding one, it returns True plus a short explanation including the index of the carbon.
If none is found (or the SMILES is invalid) the function returns False with a reason.
"""

from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    
    A primary alcohol must have at least one –OH group attached to a saturated (sp3) carbon 
    that has either:
      • three hydrogens attached (making it a CH3–OH), or
      • two hydrogens attached and exactly one carbon neighbor (making it a RCH2–OH).
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if at least one primary alcohol group is found, False otherwise.
        str: Explanation for the classification decision.
    """
    # Try to parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are reliable.
    mol = Chem.AddHs(mol)
    
    # Iterate over atoms to search for possible hydroxyl groups.
    for atom in mol.GetAtoms():
        # We are interested only in oxygen atoms.
        if atom.GetSymbol() != "O":
            continue

        # For a simple hydroxyl (–OH) group the oxygen should have exactly one heavy-atom neighbor.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue  # Skip if oxygen is attached to more than one heavy atom.
        
        # The one heavy neighbor should be a carbon.
        carbon = heavy_neighbors[0]
        if carbon.GetSymbol() != "C":
            continue

        # Ensure the carbon is saturated (sp3).
        if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # Count the number of hydrogen atoms attached to the carbon.
        num_H = sum(1 for nbr in carbon.GetNeighbors() if nbr.GetSymbol() == "H")
        # Count the number of carbon neighbors (exclude the oxygen that is attached).
        num_C_neighbors = sum(1 for nbr in carbon.GetNeighbors() if nbr.GetSymbol() == "C")
        
        # Case 1: A methanol moiety (CH3–OH) will have three H atoms.
        if num_H == 3:
            return True, f"Found primary alcohol group: CH3-OH at carbon atom index {carbon.GetIdx()}."
        
        # Case 2: A primary alcohol RCH2–OH should have exactly two hydrogens and one C neighbor.
        if num_H == 2 and num_C_neighbors == 1:
            return True, f"Found primary alcohol group: RCH2-OH at carbon atom index {carbon.GetIdx()}."
    
    return False, "No primary alcohol group found"

# For testing purposes (remove or comment out for production use):
if __name__ == '__main__':
    # Example SMILES for primary alcohols
    test_smiles = [
        "Nc1nc(=O)[nH]cc1CO",  # 5-(hydroxymethyl)cytosine (should be True)
        "CCC(CC)Nc1cc(ccc1N1C(=O)CCC1(CO)CO)C(O)=O",  # Example with –(RCH2–OH)
        "OC1C=C(Cl)C(Cl)C=C1Cl",  # Secondary structure: likely False
        "O(C([H])([H])[H])[H]"  # methanol-d1 (should be True)
    ]
    for s in test_smiles:
        result, reason = is_primary_alcohol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*40}")