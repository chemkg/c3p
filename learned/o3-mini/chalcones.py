"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones (1,3-diphenylpropenone and its substituted derivatives)
A chalcone is defined as a ketone with the core α,β‐unsaturated system attached to two aromatic rings.
This program checks for the presence of the chalcone substructure,
either as Ar–CH=CH–C(=O)–Ar or as the reverse Ar–C(=O)–CH=CH–Ar.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone (or derivative) based on its SMILES string.
    A chalcone contains an alpha,beta-unsaturated ketone linker with two aromatic groups attached.
    
    This function checks for two possible substructure patterns:
       Pattern 1: Aromatic–CH=CH–C(=O)–Aromatic (e.g., c1ccccc1C=CC(=O)c2ccccc2)
       Pattern 2: Aromatic–C(=O)–CH=CH–Aromatic (the reverse ordering)
    
    Args:
        smiles (str): SMILES representation of the molecule
        
    Returns:
        bool: True if the molecule contains a chalcone-like substructure, False otherwise.
        str: A descriptive reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the chalcone substructure.
    # Here, "a" matches any aromatic atom.
    # Pattern 1: aromatic–CH=CH–C(=O)–aromatic
    chalcone_pattern1 = Chem.MolFromSmarts("a-[CH]=[CH]-[C](=O)-a")
    # Pattern 2: aromatic–C(=O)–CH=CH–aromatic (the reverse connectivity)
    chalcone_pattern2 = Chem.MolFromSmarts("a-[C](=O)-[CH]=[CH]-a")
    
    # Check if either pattern is present in the molecule
    match1 = mol.HasSubstructMatch(chalcone_pattern1)
    match2 = mol.HasSubstructMatch(chalcone_pattern2)
    
    if match1 or match2:
        return True, "Molecule contains a chalcone-like α,β-unsaturated ketone substructure with two aromatic rings."
    else:
        return False, "Molecule does not contain the chalcone core (α,β-unsaturated ketone linked to two aromatic rings) as defined."

# Example usage:
if __name__ == "__main__":
    # testing one example: cis-chalcone, SMILES: "O=C(c1ccccc1)\\C=C/c1ccccc1"
    smiles_example = "O=C(c1ccccc1)\\C=C/c1ccccc1"
    result, reason = is_chalcones(smiles_example)
    print(f"Result: {result}\nReason: {reason}")