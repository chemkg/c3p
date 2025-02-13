"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: Secondary alpha-hydroxy ketone (acyloin)

A secondary alpha-hydroxy ketone has the structural motif R–C(=O)–CH(OH)–R′,
where the CH(OH) center is secondary (i.e. it bears exactly one hydrogen)
and is bonded to a carbonyl (C(=O)) on one side and an organic (carbon) group on the other.
This motif is typically formed by the reductive coupling of two carboxylic acid groups.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone (acyloin)
    from its SMILES string.
    
    This function searches for the acyloin motif that satisfies:
      • a ketone group (C(=O)) that is directly bonded to the acyloin center,
      • a central carbon that is secondary (exactly one hydrogen, indicated with H1)
        and has three total connections (X3) – one to the carbonyl, one to an OH group,
        and one to an organic (carbon) substituent,
      • and a second carbon substituent (organyl group) attached to the acyloin center.
    
    The corresponding SMARTS used is: "[#6]-C(=O)-[C;H1;X3](O)-[#6]"
    which helps to minimize both false negatives and false positives.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains the secondary alpha-hydroxy ketone motif, False otherwise.
        str: A message explaining the result.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for a secondary alpha-hydroxy ketone (acyloin).
    # Pattern breakdown:
    #   [#6]            -> A carbon atom as the organic substituent.
    #   -C(=O)-        -> A ketone group attached to the acyloin center.
    #   [C;H1;X3](O)    -> A carbon with exactly one hydrogen (secondary) and three connections (X3)
    #                     that bears an OH group.
    #   -[#6]          -> The other organic substituent (a carbon) attached to the acyloin center.
    pattern = Chem.MolFromSmarts("[#6]-C(=O)-[C;H1;X3](O)-[#6]")
    if pattern is None:
        return False, "Error parsing SMARTS pattern"
        
    # Check if the molecule contains the acyloin motif.
    if mol.HasSubstructMatch(pattern):
        return True, "Contains secondary alpha-hydroxy ketone motif (acyloin)"
    else:
        return False, "Does not contain the secondary alpha-hydroxy ketone (acyloin) motif"

# Example usage for testing purposes:
if __name__ == "__main__":
    test_smiles = [
        "OC(C(=O)c1ccccc1)c1ccccc1",  # benzoin: should be classified as acyloin (True)
        "OC[C@H](O)C(=O)CO",          # D-erythrulose: should be classified as acyloin (True)
        "O=C1/C(=C(/O)\\C=C\\C=C\\C)/[C@H]2[C@@](O)(C)C([C@@]1(C(OCC)C2)C)=O"  # Rezishanone C: likely False
    ]
    for smi in test_smiles:
        result, reason = is_secondary_alpha_hydroxy_ketone(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")