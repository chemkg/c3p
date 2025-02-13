"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Aromatic Primary Alcohol
Definition: Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon 
            which is itself bonded to an aromatic ring.
The approach uses a SMARTS pattern that matches a CH2OH group (sp3 carbon with two hydrogens,
one hydroxyl substituent, and one aromatic substituent).
"""

from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule contains an aromatic primary alcohol group.
    A qualifying structure has a CH2OH group (primary alcohol) where the CH2 carbon is sp3,
    has exactly two hydrogens (H2), is bonded to an -OH group (ensuring the oxygen has at least one hydrogen),
    and its one heavy-atom neighbor (other substituent) is aromatic.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule contains an aromatic primary alcohol group, False otherwise.
        str: Explanation/reason for the decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to ensure hydrogen counts are correct
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for an aromatic primary alcohol: an sp3 carbon with exactly 2 Hs,
    # one -OH substituent and one aromatic substituent.
    # The SMARTS "[C;X4;H2]([OH])([a])" matches a CH2 group with one hydroxyl group and one aromatic neighbor.
    pattern = Chem.MolFromSmarts("[C;X4;H2]([OH])([a])")
    if pattern is None:
        return False, "Error in defining the SMARTS pattern"
    
    # If a match exists, then we have an aromatic primary alcohol in the molecule.
    if mol.HasSubstructMatch(pattern):
        return True, "Found an aromatic primary alcohol group (CH2OH attached directly to an aromatic ring)"
    else:
        return False, "No aromatic primary alcohol group found in the molecule"

# For testing purposes (which you may remove when using this as a module)
if __name__ == "__main__":
    test_smiles = [
        "CCCC1=C(C(O)=O)C(CO)=CO1",           # 4-(hydroxymethyl)-2-propylfuran-3-carboxylic acid
        "CCCCn1c(CO)ccc1C=O",                 # 1-butylpyrraline
        "CC1=CC=C(CO)C(C)=C1",                # 2,4-dimethylbenzyl alcohol
        "C1=CC(=CN=C1)CO",                    # 3-pyridinemethanol
        "CCCC"                               # simple alkane (should be False)
    ]
    
    for smi in test_smiles:
        result, reason = is_aromatic_primary_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")