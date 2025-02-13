"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide
Definition: An amide of a sulfonic acid RS(=O)2NR'2.
This module attempts to verify if a molecule contains a sulfonamide functional group,
i.e. a sulfonyl group (S double-bonded to two oxygens) connected via a single bond to a nitrogen,
where that nitrogen is further bonded (one or two bonds to heavy atoms) as expected.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as an amide of a sulfonic acid, with the general structure:
      RS(=O)(=O)NR'2  (i.e. the sulfonyl group attached to a nitrogen that carries one or two substituents)

    This function first uses a SMARTS pattern to locate a candidate sulfonamide moiety.
    Then it inspects the matched nitrogen atom to make sure that aside 
    from the bonded sulfur, it is attached to one or two heavy‐atoms (C, N, O, etc.)
    as expected for a sulfonamide (this helps to avoid false positives).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a sulfonamide, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a sulfonamide group.
    # This pattern requires a tetra-coordinated sulfur (atomic num 16) with two double-bonded oxygens,
    # connected via a single bond to a nitrogen.
    sulfonamide_pattern = Chem.MolFromSmarts("[#16X4](=[OX1])(=[OX1])[NX3]")
    
    # Find all substructure matches of the pattern.
    matches = mol.GetSubstructMatches(sulfonamide_pattern)
    if not matches:
        return False, "No sulfonamide group (S(=O)(=O)N) found in the molecule"
    
    # Loop through each match to filter out hits that are not valid sulfonamides.
    # The assumed order in the pattern is: index 0 = S; index 1 and 2 = O; index 3 = N.
    valid_match_found = False
    for match in matches:
        # In our SMARTS the first atom is the sulfur and the fourth is the nitrogen.
        # (Note: this ordering is usually preserved by RDKit.)
        s_idx = match[0]
        n_idx = match[3]
        
        # Retrieve the nitrogen atom.
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Count heavy-atom neighbors (ignoring hydrogens) except the sulfur in this match.
        heavy_neighbors = [nbr for nbr in n_atom.GetNeighbors() 
                           if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != s_idx]
        # For a sulfonamide, the nitrogen is expected to have one (NH) or two (N-alkyl) heavy-atom neighbors.
        if len(heavy_neighbors) in [1, 2]:
            valid_match_found = True
            break
            
    if valid_match_found:
        return True, "Molecule contains a sulfonamide group: RS(=O)(=O)NR'2"
    else:
        return False, "No valid sulfonamide group found in the molecule"
    
# Example usage:
# Uncomment the following lines to test a few examples:
# test_smiles = [
#     "Brc1ccc(c2ccccc12)S(=O)(=O)NCc1ccccn1",  # pyrabactin, expected true
#     "CC1=CC=CC=C1S(=O)(=O)N2[C@@H]3CN4C(=O)C=CC=C4[C@H]2[C@H]([C@@H]3CO)C(=O)NCCOC",  # another example
# ]
# for sm in test_smiles:
#     result, reason = is_sulfonamide(sm)
#     print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")