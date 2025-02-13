"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: Secondary alpha-hydroxy ketone (acyloin)

A secondary alpha-hydroxy ketone has the structural motif R–C(=O)–CH(OH)–R′ 
(or the reverse R–CH(OH)–C(=O)–R′). The acyloin center must be secondary, meaning
it bears exactly one hydrogen and has exactly three non-hydrogen connections. This 
motif is typically produced by the reductive coupling of two carboxylic acids.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone (acyloin)
    based on its SMILES string.
    
    Because the acyloin motif R–C(=O)–CH(OH)–R′ can be drawn in two directions,
    we search for both orders.
    
    The function first adds explicit hydrogens (necessary for correctly matching 
    the secondary center with exactly one hydrogen) and then uses two SMARTS patterns:
    
      Pattern 1: "[#6]-C(=O)-[C;H1;X3](O)-[#6]"
         captures molecules where the carbonyl group is on the left of the acyloin center.
      
      Pattern 2: "[#6]-[C;H1;X3](O)-C(=O)-[#6]"
         captures molecules where the carbonyl group is on the right of the acyloin center.
    
    The acyloin center is defined as a carbon atom bearing exactly one hydrogen (H1)
    and three heavy-atom bonds (X3) with one substituent being an OH group.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains the secondary alpha-hydroxy ketone motif, False otherwise.
        str: A message explaining the result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the secondary carbon (with one H) is represented properly.
    mol = Chem.AddHs(mol)
    
    # Define two SMARTS for the acyloin motif.
    # Pattern 1: R-C(=O)-CH(OH)-R'
    pattern1 = Chem.MolFromSmarts("[#6]-C(=O)-[C;H1;X3](O)-[#6]")
    if pattern1 is None:
        return False, "Error parsing SMARTS pattern (pattern1)."
    
    # Pattern 2: R-CH(OH)-C(=O)-R'
    pattern2 = Chem.MolFromSmarts("[#6]-[C;H1;X3](O)-C(=O)-[#6]")
    if pattern2 is None:
        return False, "Error parsing SMARTS pattern (pattern2)."
    
    # Check if either pattern exists in the molecule.
    if mol.HasSubstructMatch(pattern1) or mol.HasSubstructMatch(pattern2):
        return True, "Contains secondary alpha-hydroxy ketone motif (acyloin)"
    else:
        return False, "Does not contain the secondary alpha-hydroxy ketone (acyloin) motif"

# Example usage for testing purposes:
if __name__ == "__main__":
    # A list of example SMILES corresponding to compounds reported to belong to this class.
    test_smiles = [
        "OC(C(=O)c1ccccc1)c1ccccc1",         # benzoin: expected True
        "OC[C@@H](O)C(=O)CO",                 # D-erythrulose: expected True
        "O=C1/C(=C(/O)\\C=C\\C=C\\C)/[C@H]2[C@@](O)(C)C([C@@]1(C(OCC)C2)C)=O"  # Rezishanone C: likely False
    ]
    for smi in test_smiles:
        result, reason = is_secondary_alpha_hydroxy_ketone(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")