"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
#!/usr/bin/env python
"""
Classifies: alpha-amino acid ester 
Definition:
  "The amino acid ester derivative obtained the formal condensation of an alpha-amino acid with an alcohol."
  
This script tests for the presence of the key substructure:
   R—C*(N)C(=O)O—R′
where the alpha carbon (C*) carries an amino group and its carboxyl group is esterified.
"""

from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    
    We expect to see a motif corresponding to an alpha-amino acid unit (R-CH(N)-COOH)
    where the carboxyl group is esterified (i.e. -COO-R) rather than free. The SMARTS
    pattern below tries to capture the fragment:
   
       [CX3]([NX3])[C](=O)O[C]
    
    This pattern matches an sp3 carbon bonded to a nitrogen and to a carbonyl carbon
    that is esterified via an oxygen to another carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an alpha-amino acid ester, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for an alpha-amino acid ester fragment.
    # This pattern looks for a carbon ([CX3]) that has a nitrogen substituent ([NX3])
    # and is attached to a carbonyl which is esterified (C(=O)O[C]).
    aa_ester_smarts = "[CX3]([NX3])[C](=O)O[C]"
    pattern = Chem.MolFromSmarts(aa_ester_smarts)
    if pattern is None:
        return False, "Failed to create SMARTS pattern"
    
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains an alpha-amino acid ester motif"
    else:
        return False, "Molecule does not contain a clear alpha-amino acid ester motif"

# For testing:
if __name__ == "__main__":
    # Example: methyl glycinate (an alpha-amino acid ester)
    test_smiles = "COC(=O)CN" 
    result, reason = is_alpha_amino_acid_ester(test_smiles)
    print(f"SMILES: {test_smiles}")
    print(f"Classification: {result}")
    print(f"Reason: {reason}")