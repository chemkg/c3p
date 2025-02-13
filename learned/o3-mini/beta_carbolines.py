"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: beta-carbolines
Definition: Any molecule containing a pyridoindole core (i.e. a beta-carboline skeleton)
           and their hydrogenated derivatives.
           
Revisions in this approach:
 – Instead of a single SMARTS we now use two refined patterns:
    • An aromatic beta-carboline pattern: "c1ccc2c(c1)[nH]c3cnccc23"
    • A tetrahydro-beta-carboline (hydrogenated) variant: "c1ccc2c(c1)[nH]C3CCncc23"
 – We require that the recognizable nitrogen appears as [nH] (an unsubstituted nitrogen) to 
   avoid spurious matches from canthinone-like structures.
 – We try each pattern; if any pattern is found the molecule is classified as in the beta-carboline class.
 
Note: This approach might not be perfect but appears to reduce both false negatives and false positives.
"""

from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule belongs to the beta-carbolines class.
    
    The classification is based on matching one of two SMARTS:
      - "c1ccc2c(c1)[nH]c3cnccc23" for an aromatic beta-carboline (pyridoindole) scaffold.
      - "c1ccc2c(c1)[nH]C3CCncc23" for a typical tetrahydro-beta-carboline variant.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str):
            True and a reason if the molecule is recognized as a beta-carboline,
            or False with a reason otherwise.
    """
    
    # Convert the SMILES string to an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS patterns for the beta-carboline core:
    patterns = [
        ("aromatic beta-carboline (pyridoindole) skeleton", "c1ccc2c(c1)[nH]c3cnccc23"),
        ("tetrahydro-beta-carboline skeleton variant", "c1ccc2c(c1)[nH]C3CCncc23")
    ]
    
    # Try each pattern
    for description, smart in patterns:
        pattern = Chem.MolFromSmarts(smart)
        if pattern is None:
            continue  # skip if SMARTS is somehow invalid
        
        if mol.HasSubstructMatch(pattern):
            return True, f"Molecule contains {description}."
    
    # If none match, then the beta-carboline substructure was not recognized.
    return False, "Molecule does not contain a recognizable beta-carboline (pyridoindole) skeleton."

# (Optional) Example usage:
if __name__ == '__main__':
    # Example: Eudistomin N (a true positive) 
    test_smiles = "Brc1ccc2[nH]c3cnccc3c2c1"
    result, reason = is_beta_carbolines(test_smiles)
    print(result, reason)