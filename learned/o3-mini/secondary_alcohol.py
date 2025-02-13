"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
Definition: A secondary alcohol is a compound in which a hydroxy group (-OH) 
is attached to a saturated (sp3) carbon atom that has exactly two carbon substituents and one hydrogen.
This function uses RDKit to parse the SMILES string, add explicit hydrogens, and search for the secondary
alcohol pattern using a SMARTS substructure.
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule contains at least one secondary alcohol group.
    
    A secondary alcohol group is defined as an sp3 (tetrahedral) carbon atom that is bonded to exactly:
      - one -OH group (where the oxygen has one hydrogen),
      - two other carbon atoms,
      - one hydrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if at least one secondary alcohol group is found, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For accurate hydrogen counts, add explicit hydrogens.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a secondary alcohol.
    # The pattern "[C;X4;H1]([C;X4])([C;X4])[O;H1]" means:
    #  • A tetrahedral carbon (C with X4) that has exactly one hydrogen (H1)
    #  • that is connected to two other sp³ carbons ([C;X4]) and
    #  • an oxygen that in turn has exactly one hydrogen ([O;H1]) (i.e. –OH).
    sec_alcohol_smarts = "[C;X4;H1]([C;X4])([C;X4])[O;H1]"
    substruct = Chem.MolFromSmarts(sec_alcohol_smarts)
    if substruct is None:
        # In an unlikely event the SMARTS did not parse.
        return False, "Failed to parse secondary alcohol SMARTS pattern."
    
    # Look for substructure matches.
    matches = mol.GetSubstructMatches(substruct)
    if matches:
        # For reporting, we report the first matching carbon (atom index returned by the match).
        # Note: the match tuple ordering is as in the SMARTS string.
        carbon_idx = matches[0][0]
        return True, f"Molecule contains a secondary alcohol group at carbon atom index {carbon_idx}."
    
    return False, "No secondary alcohol group found."

# Example usage (you can remove these lines when integrating into your pipeline):
if __name__ == "__main__":
    # Examples provided in the request:
    examples = [
        ("O[C@H]([C@H](O)CO)CCO", "2-deoxy-D-erythro-pentitol"),
        ("C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@](C)(CCC=C(C)C)O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O", "ginsenoside Re"),
        ("O[C@@H](CC\\C=C\\c1ccccc1)CCc1ccc(O)c(O)c1", "(-)-(3S)-1-(3,4-dihydroxyphenyl)-7-phenyl-(6E)-6-hepten-3-ol"),
        ("C[C@H](O)[C@H](C)O", "(S,S)-butane-2,3-diol"),
    ]
    for smi, name in examples:
        flag, reason = is_secondary_alcohol(smi)
        print(f"SMILES: {smi}\nName: {name}\nClassification: {flag}\nReason: {reason}\n")