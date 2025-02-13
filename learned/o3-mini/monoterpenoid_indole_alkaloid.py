"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: Monoterpenoid Indole Alkaloids
Definition: A terpenoid indole alkaloid which is biosynthesised from L‐tryptophan 
            and diisoprenoid (usually secolaganin) building blocks.
Heuristic:
  - The molecule must be a valid one.
  - It must contain an indole moiety. Here we check using two related SMARTS patterns.
  - After ‘removing’ the indole, the remaining fragment (ideally the diisoprenoid part)
    is expected to contain roughly 8–15 carbon atoms.
  - Overall, the total carbon count should be in the range ~20–35, at least 2 nitrogen atoms,
    and molecular weight roughly between 250 and 600 Da.
  - Additionally, we check that the molecule is sufficiently polycyclic (≥3 rings).
Note:
  This is only one heuristic approach. Many monoterpenoid indole alkaloids have very complex 
  skeletons and biosynthetic history so the classification (and performance metrics) may vary.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    
    The heuristic algorithm applies the following criteria:
      1. The SMILES parses to a valid molecule.
      2. The molecule must contain an indole substructure. We search using two related SMARTS.
      3. On “deleting” the indole substructure, the largest remaining fragment should have a 
         carbon count characteristic of a monoterpene (roughly 8 to 15 carbon atoms).
      4. In addition, the overall molecular properties are checked:
         - Total carbon count between 20 and 35.
         - At least 2 nitrogen atoms.
         - Molecular weight between 250 and 600 Da.
         - At least 3 rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True plus a message if the molecule meets the criteria, otherwise False 
                     and a message explaining why.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular properties.
    # Total carbon atom count.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (20 <= total_c <= 35):
        return False, f"Total carbon count of {total_c} outside expected range (20-35)"
    
    # Total nitrogen count.
    total_n = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if total_n < 2:
        return False, f"Found only {total_n} nitrogen atoms; expect at least 2"
    
    # Molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mw <= 600):
        return False, f"Molecular weight of {mw:.1f} Da outside expected range (250-600 Da)"
    
    # Check number of rings.
    n_rings = Chem.GetSSSR(mol)
    if n_rings < 3:
        return False, f"Found only {n_rings} rings; expected at least 3 for complex alkaloids"
    
    # Look for an indole substructure. 
    # We check two SMARTS patterns (they are very similar but may catch slightly different cases).
    indole_smarts_list = [
        "c1ccc2c(c1)[nH]c(c2)",  # typical indole
        "c1cc2c(c1)[nH]cc2"       # slight variant
    ]
    indole_found = False
    for smarts in indole_smarts_list:
        indole_pat = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(indole_pat):
            indole_found = True
            chosen_indole = indole_pat  # use this pattern for removal later
            break
    if not indole_found:
        return False, "No indole moiety found"
    
    # Delete the indole substructure from the molecule.
    # This will (heuristically) leave the remaining fragment(s) that might correspond 
    # to the diisoprenoid (monoterpene) unit.
    mol_no_indole = Chem.DeleteSubstructs(mol, chosen_indole)
    # Sometimes DeleteSubstructs leaves “dummy” atoms or unsanitized molecules; try sanitizing.
    Chem.SanitizeMol(mol_no_indole, catchErrors=True)
    
    # Decompose remaining molecule into fragments (if disconnected).
    frags = Chem.GetMolFrags(mol_no_indole, asMols=True)
    if not frags:
        return False, "Indole found but could not extract residual fragment"
    
    # Choose the largest fragment in terms of carbon count.
    max_c = 0
    for frag in frags:
        frag_c = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if frag_c > max_c:
            max_c = frag_c
    # We expect the diisoprenoid (monoterpene) fragment to have roughly 8-15 carbon atoms.
    if not (8 <= max_c <= 15):
        return False, f"Largest residual fragment has {max_c} carbons; not typical for a diisoprenoid unit"
    
    return True, ("Contains indole moiety, evidence of a diisoprenoid fragment (after indole removal), and " +
                  "typical molecular properties for a monoterpenoid indole alkaloid")

# Example usage (for testing); you can add more test cases as needed.
if __name__ == "__main__":
    test_smiles = "CC[C@]1(CN2CCC=3C4=CC(=C(C=C4NC3[C@@]2(C[C@@]1(CCO)[H])[H])OC)OC)[H]"  # Example: Ochropposinine
    result, reason = is_monoterpenoid_indole_alkaloid(test_smiles)
    print(result, reason)