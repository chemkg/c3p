"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: beta-carbolines
Definition: Any pyridoindole containing a beta-carboline skeleton (i.e. an indole fused with a pyridine)
           and hydrogenated derivatives thereof.
           
This implementation first checks for a substructure match to one of two SMARTS patterns:
  (a) an aromatic beta-carboline (pyrido[3,4-b]indole) pattern, and 
  (b) a tetrahydro-beta-carboline pattern.
If neither of these is found, it further tries to find an indole fragment fused with a pyridine fragment,
i.e. look for an indole ring ("c1ccc2[nH]cc2c1") and a pyridine ring ("c1ncccc1") that share at least two atoms.
This multi‐pronged approach attempts to catch both decorated and partially saturated variants.
"""

from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule belongs to the beta-carbolines class.
    
    Beta-carbolines (and their hydrogenated derivatives) are defined here as molecules that
    contain a pyridoindole core. In fully aromatic beta-carbolines this core is that of a
    pyrido[3,4-b]indole. However, many beta-carboline derivatives are decorated and/or partially
    saturated. To improve sensitivity while avoiding many false positives our approach is:
      (1) Use two SMARTS patterns:
             - aromatic_pattern: matching an aromatic beta-carboline skeleton.
             - tetrahydro_pattern: matching a frequently encountered tetrahydro-beta-carboline version.
      (2) If neither fires then check for the following:
             - an indole fragment ("c1ccc2[nH]cc2c1"), and 
             - a pyridine fragment ("c1ncccc1")
          that overlap (share two or more atoms). This “fused‐ring” criterion is intended to capture
          cases where the beta-carboline core is embedded in larger or spirocyclic systems.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True and an explanation if the molecule is classified as a beta-carboline
                     (or a common hydrogenated variant); otherwise False and a reason.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Attempt to match two SMARTS queries.
    # SMARTS for fully aromatic beta-carboline (pyrido[3,4-b]indole)
    aromatic_smarts = "c1ccc2c(c1)[nH]c3cnccc23"
    aromatic_pattern = Chem.MolFromSmarts(aromatic_smarts)
    
    # SMARTS for a common tetrahydro-beta-carboline variant.
    # (Here we allow for the partial saturation in the pyridine part.)
    tetrahydro_smarts = "c1ccc2c(c1)[nH]CC2"
    tetrahydro_pattern = Chem.MolFromSmarts(tetrahydro_smarts)
    
    if mol.HasSubstructMatch(aromatic_pattern):
        return True, "Molecule contains an aromatic beta-carboline (pyridoindole) skeleton."
    if mol.HasSubstructMatch(tetrahydro_pattern):
        return True, "Molecule contains a tetrahydro-beta-carboline skeleton variant."
    
    # (2) As a backup, search for an indole fragment fused with a pyridine ring.
    # Define an indole SMARTS (benzene fused to a pyrrole).
    indole_smarts = "c1ccc2[nH]cc2c1"
    indole_pattern = Chem.MolFromSmarts(indole_smarts)
    # Define a pyridine SMARTS.
    pyridine_smarts = "c1ncccc1"
    pyridine_pattern = Chem.MolFromSmarts(pyridine_smarts)
    
    indole_matches = mol.GetSubstructMatches(indole_pattern)
    pyridine_matches = mol.GetSubstructMatches(pyridine_pattern)
    
    # Check if any indole and pyridine substructure share at least 2 atoms.
    for indole_atoms in indole_matches:
        set_indole = set(indole_atoms)
        for pyridine_atoms in pyridine_matches:
            # Calculate the intersection (fused ring criterion).
            if len(set_indole.intersection(pyridine_atoms)) >= 2:
                return True, ("Molecule contains a beta-carboline (pyridoindole) core as demonstrated by "
                              "the fusion of an indole ring with a pyridine ring.")
    
    # If none of the above tests fire then we do not recognize a beta-carboline skeleton.
    return False, "Molecule does not contain a recognizable beta-carboline (pyridoindole) skeleton."

# (Optional) Example usage:
if __name__ == '__main__':
    # Example SMILES (one of the provided true positives):
    test_smiles = "Brc1ccc2[nH]c3cnccc3c2c1"  # Eudistomin N
    result, reason = is_beta_carbolines(test_smiles)
    print(result, reason)