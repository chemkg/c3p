"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: beta-carbolines (i.e. any pyridoindole containing a beta-carboline skeleton and its hydrogenated derivatives)
Note: This code uses two substructure patterns. The first (aromatic_pattern) is meant to capture the fully aromatic
beta-carboline structure (pyrido[3,4-b]indole). The second (saturated_pattern) is a simplified example to capture a frequently
found hydrogenated variant, such as tetrahydro-beta-carboline. Depending on the diversity of beta-carboline derivatives,
these patterns may need further refinement.
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule belongs to the beta-carbolines class based on the presence of a beta-carboline skeleton
    (pyridoindole) or its hydrogenated derivatives.
    
    Beta-carbolines are defined as any pyridoindole in which an indole ring (benzene fused to a pyrrole)
    is fused with an additional pyridine ring. Many natural and synthetic beta-carbolines are decorated, and some
    may have partial saturation. Here we implement two SMARTS queries:
      1. An aromatic beta-carboline query.
      2. A relaxed (partially hydrogenated) beta-carboline query.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule matches a beta-carboline skeleton or its common hydrogenated variant.
        str: Explanation of the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an aromatic beta-carboline skeleton.
    # This pattern broadly corresponds to a pyrido[3,4-b]indole.
    aromatic_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3cnccc23")
    
    # Define a SMARTS pattern for a commonly encountered hydrogenated derivative.
    # For example, "2,3,4,9-tetrahydro-beta-carboline" often has the pyridine part partially saturated.
    # This pattern is only one representative example and may not cover all cases.
    saturated_pattern = Chem.MolFromSmarts("[CH2]1CNCC2=[c]1[c]3cccc3N2")
    
    # Check if the molecule contains either the aromatic beta-carboline or the saturated (hydrogenated) derivative.
    if mol.HasSubstructMatch(aromatic_pattern):
        return True, "Molecule contains an aromatic beta-carboline (pyridoindole) skeleton."
    if mol.HasSubstructMatch(saturated_pattern):
        return True, "Molecule contains a hydrogenated beta-carboline (pyridoindole) skeleton variant."
    
    # If neither pattern is found, then the molecule does not belong to this class.
    return False, "Molecule does not contain a recognizable beta-carboline skeleton."