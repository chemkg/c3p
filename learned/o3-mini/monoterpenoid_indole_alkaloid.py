"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: Monoterpenoid Indole Alkaloids
Definition: A terpenoid indole alkaloid which is biosynthesised from L‐tryptophan 
            and diisoprenoid (usually secolaganin) building blocks.
Heuristic:
  - The molecule must be a valid one.
  - It must contain an indole moiety. Two related SMARTS patterns are used.
  - After “removing” the indole, the largest remaining fragment (if any)
    should have roughly 8–15 carbon atoms to constitute the diisoprenoid part.
  - Additionally, the overall molecular properties are checked:
       * Total carbon count between 20 and 35.
       * At least 2 nitrogen atoms.
       * Molecular weight between 250 and 600 Da.
       * Sufficient polycyclicity (≥3 rings).
Note:
  This heuristic is not definitive given the complexity of monoterpenoid indole alkaloid scaffolds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    
    The heuristic algorithm applies the following criteria:
      1. The SMILES parses to a valid molecule.
      2. The molecule must contain an indole substructure – checked with two SMARTS patterns.
      3. After “deleting” the indole substructure, the largest remaining fragment should be 
         consistent with a diisoprenoid (monoterpene) unit (roughly 8-15 carbon atoms).
      4. In addition, overall molecular properties are checked:
         - Total carbon count between 20 and 35.
         - At least 2 nitrogen atoms.
         - Molecular weight between 250 and 600 Da.
         - At least 3 rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple where the first value is True if the molecule meets the criteria,
                     otherwise False, and the second value gives the reasoning.
    """
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure a clean molecule by removing explicit hydrogens.
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
    
    # Count rings using rdMolDescriptors function which reliably returns an integer.
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3:
        return False, f"Found only {n_rings} rings; expected at least 3 for complex alkaloids"
    
    # Search for an indole substructure using two similar SMARTS patterns.
    indole_smarts_list = [
        "c1ccc2c(c1)[nH]c(c2)",  # standard indole
        "c1cc2c(c1)[nH]cc2"       # slight variant
    ]
    indole_found = False
    chosen_indole = None
    for smarts in indole_smarts_list:
        indole_pat = Chem.MolFromSmarts(smarts)
        if indole_pat is None:
            continue  # skip if pattern fails to generate
        if mol.HasSubstructMatch(indole_pat):
            indole_found = True
            chosen_indole = indole_pat  # remember this pattern for deletion
            break
    if not indole_found:
        return False, "No indole moiety found"
    
    # Remove the indole substructure.
    mol_no_indole = Chem.DeleteSubstructs(mol, chosen_indole)
    # Remove any dummy atoms and explicit hydrogens left over and sanitize molecule.
    mol_no_indole = Chem.RemoveHs(mol_no_indole)
    try:
        Chem.SanitizeMol(mol_no_indole, catchErrors=True)
    except Exception as e:
        # If sanitize fails, return an appropriate error.
        return False, f"Sanitization after indole removal failed: {e}"
    
    # Decompose molecule into fragments.
    frags = Chem.GetMolFrags(mol_no_indole, asMols=True)
    if not frags:
        return False, "Indole found but no residual fragment could be extracted"
    
    # Choose the largest fragment by number of carbon atoms.
    max_c = 0
    for frag in frags:
        frag_c = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if frag_c > max_c:
            max_c = frag_c
    # Expect the diisoprenoid (monoterpene) fragment to have roughly 8-15 carbon atoms.
    if not (8 <= max_c <= 15):
        return False, f"Largest residual fragment has {max_c} carbons; not typical for a diisoprenoid unit"
    
    return True, ("Contains indole moiety, evidence of a diisoprenoid fragment (after indole removal), and " +
                  "typical molecular properties for a monoterpenoid indole alkaloid")


# Example usage (for testing)
if __name__ == "__main__":
    # Test with Ochropposinine as an example.
    test_smiles = "CC[C@]1(CN2CCC=3C4=CC(=C(C=C4NC3[C@@]2(C[C@@]1(CCO)[H])[H])OC)OC)[H]"
    result, reason = is_monoterpenoid_indole_alkaloid(test_smiles)
    print(result, reason)