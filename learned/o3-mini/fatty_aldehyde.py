"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehydes
Definition: “An aldehyde formally arising from reduction of the carboxylic acid group 
of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.”
This improved routine uses additional heuristics:
  - The molecule must be successfully parsed.
  - It must contain only atoms typical for a fatty aldehyde (C, H, O) and no rings.
  - It must have a minimum number of carbons and a minimum molecular weight.
  - There must be at least one aldehyde group ([CX3H1](=O)),
    and that aldehyde must be terminal – meaning that its carbonyl carbon is attached 
    to exactly one carbon atom. Further, if we “walk” from that neighbor along carbon–carbon 
    bonds (ignoring double/triple bond differences) we require a linear segment of a certain length.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is defined as arising from a fatty acid reduction and
    should feature a terminal aldehyde group attached to a long, unbranched (or nearly unbranched)
    aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation of the decision.
    """
    # Try parsing the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Heuristic 1: Molecule should contain only C, H, and O.
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, ("Contains atoms other than C, H, and O "
                           "(e.g. heteroatoms, metals or others) which is not typical for a fatty aldehyde")

    # Heuristic 2: Fatty aldehydes originate from fatty acids, typically acyclic. 
    # Reject molecules with rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings which is not typical for a fatty aldehyde"
    
    # Heuristic 3: Check that the molecule has a minimum number of carbon atoms 
    # and a minimum molecular weight (exclude small aldehydes not derived from fatty acids).
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, f"Not enough carbon atoms ({c_count}) for a fatty aldehyde"
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a fatty aldehyde"
    
    # Heuristic 4: Look for an aldehyde group with the SMARTS pattern [CX3H1](=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not matches:
        return False, "No aldehyde group found"

    # We require that at least one aldehyde group be terminal.
    # For a terminal aldehyde, the carbonyl carbon should have exactly one neighboring carbon.
    # Then, starting from that neighbor, we “walk” along the chain (only through carbons)
    # and require a minimal linear chain length.
    # (This helps to discriminate simple aldehydes and complex molecules with extraneous groups.)
    
    # Set the minimum contiguous carbon chain length (starting from the neighbor of the aldehyde)
    # that we will accept as representing a 'fatty' chain.
    MIN_CHAIN_LENGTH = 3  # adjust as needed (e.g. nonanal has an 8-carbon chain overall, but even 3 in a contiguous chain may be acceptable if MW filter is met)
    
    terminal_aldehyde_found = False
    for match in matches:
        # match[0] is the carbonyl carbon in the aldehyde group.
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        # Count neighboring carbons (ignoring oxygen neighbors)
        neighbor_carbons = [nbr for nbr in aldehyde_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(neighbor_carbons) != 1:
            # Not terminal if the carbonyl C has more than one carbon neighbor.
            continue
        # Starting from the only carbon neighbor, try to count a contiguous linear chain.
        chain_length = 1
        prev_atom = aldehyde_c
        current = neighbor_carbons[0]
        # We iterate along the chain while there is exactly one continuation.
        # (If a branching occurs, this walk stops.)
        while True:
            # Find carbon neighbors of 'current', excluding the one we came from.
            cont_neighbors = [nbr for nbr in current.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx()]
            if len(cont_neighbors) == 1:
                chain_length += 1
                prev_atom = current
                current = cont_neighbors[0]
            else:
                break
        # Check if this chain length meets our threshold.
        if chain_length >= MIN_CHAIN_LENGTH:
            terminal_aldehyde_found = True
            break

    if not terminal_aldehyde_found:
        return (False, "Aldehyde group is present but does not appear terminal "
                       "or is attached to an insufficiently long unbranched aliphatic chain")
    
    return True, "Contains a terminal aldehyde group attached to a long, unbranched aliphatic chain and qualifies as a fatty aldehyde"

# Example usage (uncomment for testing):
# test_smiles = "O=CCCCCCCCCC/C=C/CC"  # 11E-Tetradecenal
# print(is_fatty_aldehyde(test_smiles))