"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: Secondary alpha-hydroxy ketone (acyloin)

A secondary alpha-hydroxy ketone has the structural motif R–C(=O)–CH(OH)–R′,
where the CH(OH) carbon (the alpha carbon) is secondary (i.e. it has exactly one hydrogen)
and is connected to a carbonyl group on one side and an organic (carbon) substituent on the other.
This motif is typically generated from the reductive coupling of two carboxylic acid groups.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone (acyloin) based on its SMILES string.
    
    It searches for the substructure fragment corresponding to C(=O)[CH1]([OX2H])[#6],
    which represents a ketone group adjacent to a CH(OH) unit where the hydroxyl-bearing
    carbon has exactly one hydrogen and is linked to another carbon (organyl group).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains the secondary alpha-hydroxy ketone motif, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to generate a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for secondary alpha-hydroxy ketone:
    # C(=O)[CH1]([OX2H])[#6]
    # Breakdown:
    #   C(=O)      => Ketone carbonyl group.
    #   [CH1]      => A carbon with exactly one hydrogen.
    #   ([OX2H])   => An –OH group (oxygen with two connections and an H).
    #   [#6]       => An attached carbon (an organic substituent).
    pattern = Chem.MolFromSmarts("C(=O)[CH1]([OX2H])[#6]")
    if pattern is None:
        return False, "Error parsing SMARTS pattern"
    
    # Perform the substructure search.
    if mol.HasSubstructMatch(pattern):
        return True, "Contains secondary alpha-hydroxy ketone motif (acyloin)"
    else:
        return False, "Does not contain the secondary alpha-hydroxy ketone (acyloin) motif"
        
# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        "OC(C(=O)c1ccccc1)c1ccccc1",  # benzoin - should return True
        "OC[C@H](O)C(=O)CO"          # D-erythrulose - may not match the strict ketone pattern
    ]
    for smi in test_smiles:
        result, reason = is_secondary_alpha_hydroxy_ketone(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")