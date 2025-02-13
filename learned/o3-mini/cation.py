"""
Classifies: CHEBI:36916 cation
"""
#!/usr/bin/env python
"""
Classifies: cation
Definition:
  A monoatomic or polyatomic species having one or more elementary charges of the proton.
  This function uses the overall net formal charge of the molecule; if it is positive (greater than 0),
  the molecule is classified as a cation.
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation must carry an overall net positive charge.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a cation, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the overall net formal charge by summing the formal charges of all the atoms
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    # Check if the molecule has a positive net charge (i.e. one or more elementary charges of the proton)
    if net_charge > 0:
        return True, f"Molecule has a net positive charge of {net_charge}."
    else:
        return False, f"Molecule has a net charge of {net_charge}, thus it is not a cation."
        
# Example usage (you can remove or comment out these lines when using this module as an imported function)
if __name__ == "__main__":
    test_smiles_list = [
        "[Co+]",                  # Cobalt(1+) should be classified as a cation
        "C(CCC)CC[NH3+]",         # hexan-1-aminium should be classified as a cation
        "CC(=O)O",                # Acetic acid (neutral) should not be classified as a cation
    ]
    
    for smiles in test_smiles_list:
        result, reason = is_cation(smiles)
        print(f"SMILES: {smiles} -> Cation: {result} ({reason})")