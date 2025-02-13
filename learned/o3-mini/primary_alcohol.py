"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: Primary Alcohol
Definition: A primary alcohol is a compound in which a hydroxy group (-OH) 
is attached to a saturated carbon atom which has either three hydrogen atoms 
attached to it (e.g., CH3OH) or only one other carbon atom and two hydrogen atoms 
attached to it (e.g., RCH2OH).

The function is_primary_alcohol takes a SMILES string as input and returns a 
boolean along with a reason for the classification.
"""

from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    
    A primary alcohol must have at least one -OH group attached to a saturated 
    (sp3) carbon that has either:
      - three hydrogens attached (making it a CH3 group) or 
      - two hydrogens attached and exactly one carbon neighbor (making it an RCH2 group).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule has at least one primary alcohol group, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms to look for hydroxyl groups.
    for atom in mol.GetAtoms():
        # Look for oxygen atoms that could form an -OH group.
        if atom.GetSymbol() != "O":
            continue
        
        # Check if this oxygen is part of an -OH group.
        # In RDKit, implicit hydrogens are not listed as neighbors, so we consider the heavy atom neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue  # Not a simple -OH (could be part of an ether or deprotonated species)
        
        # Get the single heavy neighbor; this should be the carbon to which the -OH is attached.
        carbon = heavy_neighbors[0]
        if carbon.GetSymbol() != "C":
            continue  # Not attached to carbon; we need an alcohol group
        
        # Check that the carbon is sp3 (saturated) 
        if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Check the total number of hydrogens attached to this carbon (implicit+explicit).
        num_hydrogens = carbon.GetTotalNumHs()
        # Count how many neighbors (heavy atoms) are carbons (excluding the oxygen we're considering)
        num_c_neighbors = sum(1 for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 6)
        
        # There are two valid cases for primary alcohol:
        # Case 1: Carbon has three hydrogens (CH3-OH); note that then the only heavy neighbor is the O.
        # Case 2: Carbon has two hydrogens and exactly one carbon neighbor (RCH2-OH).
        if num_hydrogens == 3:
            return True, f"Found primary alcohol group: CH3-OH at carbon atom index {carbon.GetIdx()}."
        if num_hydrogens == 2 and num_c_neighbors == 1:
            return True, f"Found primary alcohol group: RCH2-OH at carbon atom index {carbon.GetIdx()}."
    
    return False, "No primary alcohol group found"