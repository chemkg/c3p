"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: Unsaturated Fatty Acid
Definition: A fatty acid (i.e. a molecule with exactly one terminal carboxyl 
group on an aliphatic (acyclic) chain) that has at least one C=C or C#C bond.
Note: This heuristic rejects molecules with multiple acid groups or acid groups 
embedded inside rings.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    A valid unsaturated fatty acid must (a) contain exactly one terminal (acyclic) carboxyl group, 
    and (b) have at least one carbon–carbon double or triple bond.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an unsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that covers both neutral and deprotonated acid groups.
    # The pattern returns a match with the acid carbon as the first atom.
    acid_smarts = "[$(C(=O)O),$(C(=O)[O-])]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Filter matches: the acid carbon (match[0]) should be acyclic and attached to exactly one carbon
    terminal_acid_matches = []
    for match in acid_matches:
        acid_c = mol.GetAtomWithIdx(match[0])
        # Reject if the acid carbon is in a ring (we require a terminal acid)
        if acid_c.IsInRing():
            continue
        # Count carbon neighbors of the acid carbon
        n_carbon_neighbors = sum(1 for nb in acid_c.GetNeighbors() if nb.GetAtomicNum() == 6)
        if n_carbon_neighbors == 1:
            terminal_acid_matches.append(match)
    
    if len(terminal_acid_matches) != 1:
        return False, ("Molecule does not have exactly one terminal carboxylic acid group "
                       f"(found {len(terminal_acid_matches)}), so not a fatty acid")
    
    # Check for unsaturation: at least one carbon-carbon double or triple bond.
    # Using SMARTS patterns for C=C and C#C bonds.
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    triple_bond_pattern = Chem.MolFromSmarts("[#6]#[#6]")
    
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    has_triple_bond = mol.HasSubstructMatch(triple_bond_pattern)
    
    if has_double_bond or has_triple_bond:
        # Prepare a descriptive message for the type of unsaturation found.
        if has_double_bond and has_triple_bond:
            unsat_detail = "carbon-carbon double and triple bond(s)"
        elif has_double_bond:
            unsat_detail = "carbon-carbon double bond(s)"
        else:
            unsat_detail = "carbon-carbon triple bond(s)"
        return True, f"Molecule is a fatty acid and contains unsaturation: {unsat_detail}"
    else:
        return False, ("Molecule appears to be a fatty acid (has a terminal carboxyl group) "
                       "but does not contain any C=C or C#C bonds (unsaturation)")
                       
# Example usage:
if __name__ == "__main__":
    # test with one of the provided examples: trans-vaccenic acid
    test_smiles = "CCCCCC\\C=C\\CCCCCCCCCC(O)=O"
    result, reason = is_unsaturated_fatty_acid(test_smiles)
    print(result, reason)