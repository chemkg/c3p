"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehydes
Definition: “An aldehyde formally arising from reduction of the carboxylic acid group of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.”
This routine uses a series of heuristics:
  - The molecule must be successfully parsed.
  - It must have a sufficient number of carbons and a minimum molecular weight.
  - It must contain at least one aldehyde group matching the pattern [CX3H1](=O).
  - At least one aldehyde group must be terminal; that is, its carbonyl carbon is connected to only one carbon and that neighbor, when discounting the aldehyde carbon, appears to be at a chain terminus.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    It uses heuristics such as the presence of a terminal aldehyde group,
    sufficient carbon count and molecular weight, and an unbranched aliphatic chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty aldehyde, False otherwise.
        str: Explanation of the reasoning behind the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, f"Not enough carbon atoms ({c_count}) for a fatty aldehyde"

    # Check that the molecular weight is in a typical range for fatty aldehydes.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a fatty aldehyde"

    # Look for an aldehyde group. We use the SMARTS pattern for an aldehyde carbon: [CX3H1](=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not matches:
        return False, "No aldehyde group found"

    # Check if at least one aldehyde group is terminal.
    # For a terminal aldehyde, the aldehyde (carbonyl) carbon should have only one carbon neighbor.
    terminal_aldehyde_found = False
    for match in matches:
        # match[0] corresponds to the aldehyde carbon from the SMARTS
        aldehyde_atom = mol.GetAtomWithIdx(match[0])
        # Get neighboring atoms that are carbons.
        neighbor_carbons = [nbr for nbr in aldehyde_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(neighbor_carbons) != 1:
            # Not the expected connectivity for a terminal group
            continue
        # To be more confident that the aldehyde is at the end of a long chain,
        # the neighbor carbon – once we conceptually remove the aldehyde carbon – should have low connectivity (i.e. be at the end of a chain)
        neighbor = neighbor_carbons[0]
        # Get the neighbor carbons of the neighbor (excluding the aldehyde carbon)
        neighbor_connections = [nbr for nbr in neighbor.GetNeighbors() 
                                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != aldehyde_atom.GetIdx()]
        if len(neighbor_connections) <= 1:
            terminal_aldehyde_found = True
            break

    if not terminal_aldehyde_found:
        return False, "Aldehyde group is present but does not appear terminal (i.e., not at one end of a long aliphatic chain)"
    
    return True, "Contains a terminal aldehyde group attached to a sufficiently long aliphatic chain and qualifies as a fatty aldehyde"