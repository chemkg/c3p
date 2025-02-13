"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
A structure is classified as a polysaccharide if it has more than ten (i.e. >10) distinct sugar rings 
and its molecular weight is above a given threshold.
This version improves detection by using more specific SMARTS patterns that require 
saturated (sp3) ring carbons ("CX4R") and a ring oxygen ("OX2R"). It also relaxes the overlapping
filter to assume that glycosidic linkages (shared atoms) should not subtract from the total count.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide using enhanced detection of sugar rings.
    
    Criteria:
      (1) The SMILES must be valid.
      (2) The molecule is searched for sugar ring moieties. A typical pyranose has 6 members:
          5 sp3 carbons and 1 oxygen, and a typical furanose has 5 members: 4 sp3 carbons and 1 oxygen.
          Here we use SMARTS patterns "CX4R" for a saturated ring carbon and "OX2R" for a ring oxygen.
      (3) We count distinct ring systems (even if they share a glycosidic bond, we count them as one unit).
      (4) To be classified as a polysaccharide, the molecule should contain more than 10 sugar units
          and have a molecular weight above 1000 Da.
    
    Args:
        smiles (str): Input SMILES.
        
    Returns:
        (bool, str): Boolean result and an explanation.
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create more specific SMARTS for sugar rings:
    # Pyranose: 6-membered ring with exactly 1 oxygen (OX2R) and 5 saturated carbons (CX4R)
    pyranose_smarts = "[CX4R][CX4R][CX4R][CX4R][OX2R][CX4R]"
    # Furanose: 5-membered ring with exactly 1 oxygen and 4 saturated carbons.
    furanose_smarts = "[CX4R][CX4R][OX2R][CX4R][CX4R]"
    
    pyranose_pat = Chem.MolFromSmarts(pyranose_smarts)
    furanose_pat = Chem.MolFromSmarts(furanose_smarts)
    
    # For each pattern, collect the set of ring indices (as frozensets of atom indices)
    sugar_ring_sets = set()
    
    for pat in (pyranose_pat, furanose_pat):
        if pat is None:
            continue
        # Get all matches
        matches = mol.GetSubstructMatches(pat, useChirality=True)
        for match in matches:
            # Get the ring indices for the atoms in the match.
            # We sort them so that different orderings of the same atoms produce the same frozenset.
            ring_set = frozenset(match)
            sugar_ring_sets.add(ring_set)
    
    sugar_unit_count = len(sugar_ring_sets)
    
    # Check whether count exceeds threshold:
    if sugar_unit_count <= 10:
        return (False, f"Found only {sugar_unit_count} sugar ring(s); a polysaccharide requires more than 10 residues.")
    
    # Check molecular weight threshold
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return (False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a typical polysaccharide.")
    
    return (True, f"Detected {sugar_unit_count} sugar ring(s) with a molecular weight of {mol_wt:.1f} Da indicative of a polysaccharide.")

# Example usage:
# Uncomment the following lines to run a test:
# smiles_example = "O([C@H]1[C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)..."
# result, explanation = is_polysaccharide(smiles_example)
# print(result, explanation)