"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: Monoterpenoid Indole Alkaloids
Definition: A terpenoid indole alkaloid which is biosynthesised from L-tryptophan 
            and diisoprenoid (usually secolaganin) building blocks.
Heuristic:
  - The molecule must have a valid structure.
  - It must contain an indole moiety (as a proxy for the L-tryptophan‐derived part).
  - It must show evidence of a diisoprenoid (monoterpene) fragment. Here we search for a 
    typical C=C-containing fragment (e.g. a “C/C=C” motif) which often appears in such compounds.
  - Additional cutoffs (total C count, nitrogen atoms, and molecular weight) are applied.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    
    This heuristic algorithm applies the following criteria:
      1. Valid molecule.
      2. Contains an indole substructure, defined here by the SMARTS "c1ccc2c(c1)[nH]c(c2)".
      3. Contains a fragment suggestive of a diisoprenoid (terpenoid) unit. We approximate 
         this by looking for a "C/C=C" pattern.
      4. Overall molecular properties: total carbon count between 15 and 30, at least 2 nitrogen atoms,
         and a molecular weight roughly between 250 and 600 Da.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True plus a success message if the molecule meets criteria, otherwise False 
                     and a reason explaining why it failed.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for indole moiety.
    # This SMARTS should match a typical indole ring: benzene fused to a pyrrole.
    indole_smarts = "c1ccc2c(c1)[nH]c(c2)"
    indole_pat = Chem.MolFromSmarts(indole_smarts)
    if not mol.HasSubstructMatch(indole_pat):
        return False, "No indole moiety found"
    
    # Check for evidence of a terpenoid (diisoprenoid) fragment.
    # A common pattern in many examples is the C/C=C motif.
    terpenoid_smarts = "C/C=C"
    terpenoid_pat = Chem.MolFromSmarts(terpenoid_smarts)
    if not mol.HasSubstructMatch(terpenoid_pat):
        return False, "No diisoprenoid (monoterpene) fragment found"
    
    # Count total carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (15 <= c_count <= 30):
        return False, f"Carbon count of {c_count} outside typical range (15-30) for monoterpenoid indole alkaloids"
    
    # Count total nitrogen atoms.
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, f"Found only {n_count} nitrogen atoms; expect at least 2 for indole alkaloids"
    
    # Check molecular weight (rough check); many such alkaloids are between 250 and 600 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mol_wt <= 600):
        return False, f"Molecular weight of {mol_wt:.1f} Da outside expected range (250-600 Da)"
    
    return True, "Contains indole moiety, diisoprenoid fragment, and typical molecular properties for a monoterpenoid indole alkaloid"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "CC[C@]1(CN2CCC=3C4=CC(=C(C=C4NC3[C@@]2(C[C@@]1(CCO)[H])[H])OC)OC)[H]"  # Ochropposinine example
    result, reason = is_monoterpenoid_indole_alkaloid(test_smiles)
    print(result, reason)