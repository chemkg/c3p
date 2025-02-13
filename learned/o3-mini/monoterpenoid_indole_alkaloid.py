"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: Monoterpenoid Indole Alkaloids
Definition: A terpenoid indole alkaloid which is biosynthesised from L‐tryptophan 
            and diisoprenoid (usually secolaganin) building blocks.
Heuristic:
  - The molecule must be valid.
  - Overall properties: total carbon count between 20 and 35, at least 2 nitrogen atoms,
    molecular weight between 250 and 600 Da, and at least 3 rings.
  - Must contain an indole moiety. Two related SMARTS patterns are tried.
  - After removal of the indole substructure, the largest remaining fragment (by number
    of carbons) should have roughly 8–15 carbons (the diisoprenoid unit).
Note:
  Due to the inherent complexity of these structures, the heuristic includes a manual
  clean-up step (Kekulization with clearAromaticFlags) after the indole deletion to address
  aromaticity issues that led to errors in previous attempts.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    
    The heuristic algorithm applies the following criteria:
      1. The SMILES is parsed to a valid molecule which is sanitized.
      2. Overall properties are checked:
         - Total carbon count in range 20-35.
         - At least 2 nitrogen atoms.
         - Molecular weight between 250 and 600 Da.
         - At least 3 rings.
      3. The molecule must contain an indole substructure (two different SMARTS patterns are used).
      4. The indole moiety is "deleted" from the molecule and the largest remaining fragment (by carbon count)
         should have roughly 8-15 carbon atoms (indicative of the diisoprenoid part).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple where the first element is True if the molecule meets the criteria,
                     and False otherwise, with the second element providing the reasoning.
    """
    
    # Parse and sanitize the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        # Do an initial sanitization to set up aromaticity etc.
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Initial sanitization failed: {e}"
    
    # Remove explicit hydrogens to ease further processing.
    mol = Chem.RemoveHs(mol)
    
    # Check overall molecular properties.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (20 <= total_c <= 35):
        return False, f"Total carbon count of {total_c} outside expected range (20-35)"
    
    total_n = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if total_n < 2:
        return False, f"Found only {total_n} nitrogen atoms; expected at least 2"
    
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mw <= 600):
        return False, f"Molecular weight of {mw:.1f} Da outside expected range (250-600 Da)"
    
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3:
        return False, f"Found only {n_rings} rings; expected at least 3 for complex alkaloids"
    
    # Look for an indole substructure using two similar SMARTS patterns.
    indole_smarts_list = [
        "c1ccc2c(c1)[nH]c(c2)",  # standard indole
        "c1cc2c(c1)[nH]cc2"       # slight variant
    ]
    
    indole_found = False
    chosen_indole = None
    for smarts in indole_smarts_list:
        indole_pat = Chem.MolFromSmarts(smarts)
        if indole_pat is None:
            continue  # skip if pattern generation fails
        if mol.HasSubstructMatch(indole_pat):
            indole_found = True
            chosen_indole = indole_pat  # keep the pattern for later removal
            break
    if not indole_found:
        return False, "No indole moiety found"
    
    # Remove the indole substructure from the molecule.
    try:
        mol_no_indole = Chem.DeleteSubstructs(mol, chosen_indole)
        mol_no_indole = Chem.RemoveHs(mol_no_indole)
    except Exception as e:
        return False, f"Error deleting indole substructure: {e}"
    
    # Re-sanitize and attempt Kekulization to clear aromaticity inconsistencies.
    try:
        Chem.SanitizeMol(mol_no_indole, catchErrors=True)
        Chem.Kekulize(mol_no_indole, clearAromaticFlags=True)
    except Exception as e:
        return False, f"Post-deletion cleanup (sanitization/Kekulization) failed: {e}"
    
    # Decompose molecule into fragments (if any remain) after deletion.
    frags = Chem.GetMolFrags(mol_no_indole, asMols=True)
    if not frags:
        return False, "Indole found but no residual fragment could be extracted after deletion"
    
    # Choose the largest fragment by counting carbon atoms.
    max_c = 0
    for frag in frags:
        frag_c = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if frag_c > max_c:
            max_c = frag_c
            
    # Expect the diisoprenoid (monoterpene) fragment to have roughly 8-15 carbon atoms.
    if not (8 <= max_c <= 15):
        return False, f"Largest residual fragment has {max_c} carbons; not typical for a diisoprenoid unit"
    
    return True, ("Contains an indole moiety, and after its removal the remaining fragment has a carbon count " +
                  f"of {max_c}, consistent with a diisoprenoid unit, plus overall molecular properties match.")

# Example usage (for testing)
if __name__ == "__main__":
    # Example test using a provided SMILES (this one from the error report).
    test_smiles = "CCC1=C(CC(=O)OC)C(C=2C3=C(C=CN2)C4=CC=CC=C4N3)(OC)OC1=O"
    result, reason = is_monoterpenoid_indole_alkaloid(test_smiles)
    print(result, reason)