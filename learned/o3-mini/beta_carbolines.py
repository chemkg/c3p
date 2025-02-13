"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: beta-carbolines
Definition: Any molecule containing a pyridoindole core (i.e. a beta-carboline skeleton)
           and its hydrogenated (e.g. tetrahydro-beta-carboline) variants.
           
Revisions from the previous approach:
 - We now use refined SMARTS for the fused ring beta-carboline systems.
 - The aromatic pattern "n1c2ccccc2c3ccncc13" catches typical beta-carbolines irrespective of N-substitution.
 - The tetrahydro pattern "n1c2ccccc2c3CCncc13" captures common hydrogenated variants.
 - We remove the fallback that looked for overlapping indole and pyridine fragments, since that
   led to many false positives.
"""

from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule belongs to the beta-carbolines class.
    
    The classification is based on whether the molecule contains either:
      - an aromatic beta-carboline (pyridoindole) core or
      - a typical tetrahydro-beta-carboline variant.
      
    Both SMARTS patterns have been refined to capture the fused bicyclic system typical
    of beta-carbolines (9H-pyrido[3,4-b]indole) even if the lone indole nitrogen is substituted.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): 
            True and a reason if the molecule is recognized as a beta-carboline, 
            or False with a reason otherwise.
    """
    
    # Convert the SMILES to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for an aromatic beta-carboline (pyridoindole) core.
    # This pattern expects a fused bicyclic system where an indole-like ring is fused to a pyridine ring.
    aromatic_smarts = "n1c2ccccc2c3ccncc13"
    aromatic_pattern = Chem.MolFromSmarts(aromatic_smarts)
    
    # SMARTS for a typical tetrahydro-beta-carboline variant.
    # Here the pyridine part is partially saturated (the two 'C's indicate reduction).
    tetrahydro_smarts = "n1c2ccccc2c3CCncc13"
    tetrahydro_pattern = Chem.MolFromSmarts(tetrahydro_smarts)
    
    # Check if the molecule matches the aromatic beta-carboline pattern.
    if mol.HasSubstructMatch(aromatic_pattern):
        return True, "Molecule contains an aromatic beta-carboline (pyridoindole) skeleton."
    
    # Check if the molecule matches the tetrahydro-beta-carboline pattern.
    if mol.HasSubstructMatch(tetrahydro_pattern):
        return True, "Molecule contains a tetrahydro-beta-carboline skeleton variant."
    
    # If neither pattern is found, then classify as not being in the beta-carbolines class.
    return False, "Molecule does not contain a recognizable beta-carboline (pyridoindole) skeleton."

# (Optional) Example usage when the script is run directly:
if __name__ == '__main__':
    # Example SMILES from provided true positives. (Eudistomin N)
    test_smiles = "Brc1ccc2[nH]c3cnccc3c2c1"
    result, reason = is_beta_carbolines(test_smiles)
    print(result, reason)