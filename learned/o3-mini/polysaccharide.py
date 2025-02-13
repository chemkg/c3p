"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
Only structures that contain more than ten monosaccharide residues (i.e. >10 sugar rings) and a sufficiently high molecular weight are classified.
This improved version uses SMARTS patterns to detect typical sugar ring scaffolds (both pyranoses and furanoses) and counts non-overlapping matches.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    
    The detection is based on the following criteria:
      - The SMILES string must be valid.
      - The molecule is searched for sugar ring moieties. For our purposes, we assume that
        a monosaccharide ring is either a five-membered (furanose) or six-membered (pyranose) ring
        that is aliphatic, contains exactly one oxygen, and the remaining ring atoms are carbons.
        To improve detection, we search with SMARTS patterns for these motifs.
      - Only nonoverlapping matches (i.e. distinct sugar units) are counted.
      - To be classified as a polysaccharide, the molecule must have >10 sugar units and possess a 
        molecular weight above a threshold (here set to 1000 Da).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polysaccharide, False otherwise.
        str: Explanation of the classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for typical sugar rings
    # Pyranose: 6-membered ring with 1 O and 5 C (nonaromatic, fully saturated)
    pyranose_smarts = "[C;R;!$(C=*)][C;R;!$(C=*)][C;R;!$(C=*)][C;R;!$(C=*)][O;R][C;R;!$(C=*)]"
    # Furanose: 5-membered ring with 1 O and 4 C
    furanose_smarts = "[C;R;!$(C=*)][C;R;!$(C=*)][O;R][C;R;!$(C=*)][C;R;!$(C=*)]"
    
    pyranose_pat = Chem.MolFromSmarts(pyranose_smarts)
    furanose_pat = Chem.MolFromSmarts(furanose_smarts)
    
    sugar_matches = []
    
    # Get matches for pyranose
    if pyranose_pat is not None:
        matches = mol.GetSubstructMatches(pyranose_pat, useChirality=True)
        sugar_matches.extend(matches)
    
    # Get matches for furanose
    if furanose_pat is not None:
        matches = mol.GetSubstructMatches(furanose_pat, useChirality=True)
        sugar_matches.extend(matches)
    
    # To avoid counting overlapping sugar rings, we filter the matches so that each atom is used only once.
    # We use a greedy algorithm: sort matches by size (all same length here) and take matches that do not share atoms.
    used_atoms = set()
    distinct_sugar_units = []
    for match in sorted(sugar_matches, key=lambda x: len(x)):
        match_set = set(match)
        if match_set & used_atoms:
            continue
        distinct_sugar_units.append(match)
        used_atoms |= match_set
    
    sugar_unit_count = len(distinct_sugar_units)
    
    # Check sugar residue count threshold (>10)
    if sugar_unit_count <= 10:
        return False, f"Found only {sugar_unit_count} sugar ring(s); a polysaccharide requires more than 10 residues."
    
    # Check molecular weight (must be high enough for a biomacromolecule)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a typical polysaccharide."
    
    return True, f"Detected {sugar_unit_count} sugar ring(s) with a molecular weight of {mol_wt:.1f} Da indicative of a polysaccharide."

# Example usage (uncomment for testing):
# result, reason = is_polysaccharide("O(C1C(O)C(O)C(O)C1O)...")  # Provide a valid polysaccharide SMILES here
# print(result, reason)