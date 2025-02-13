"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: nonclassic icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
other than the classic icosanoids (the leukotrienes and the prostanoids).

Heuristic rules used in this improved version:
  1. Molecule must be valid.
  2. Molecule must contain only C, H, and O (icosanoids are fatty acid derivatives).
  3. Total number of carbon atoms is roughly in the range 18–26.
  4. Molecular weight should be in a window that fits eicosanoid derivatives (250–500 Da).
  5. Contains at least one carboxyl functionality – either a free acid (C(=O)[OH] or ionic form),
     or an ester (C(=O)O[C]), to allow for derivatives such as methyl esters.
  6. Contains enough oxygen atoms (>= 3).
  7. Excludes molecules containing a 5‐membered ring (marker for classic prostanoids).
  
Note: This is a heuristic approach and may not cover all edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    
    Heuristic criteria:
      - Valid molecule.
      - Contains only carbon, hydrogen, and oxygen atoms.
      - Total carbon count is between 18 and 26.
      - Molecular weight is between 250 and 500 Da.
      - Contains at least one carboxyl functionality (free acid, carboxylate, or ester) 
        as determined by SMARTS.
      - Has at least 3 oxygen atoms.
      - Does not contain any 5-membered rings.
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a nonclassic icosanoid, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that molecule contains only C, H, and O.
    # (Hydrogen atoms are implicit so we check only heavy atoms.)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 8]:
            return False, f"Contains element with atomic number {atomic_num} not allowed for icosanoids (only C, H, O permitted)"
    
    # Count total carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)
    if num_carbons < 18 or num_carbons > 26:
        return False, f"Carbon count {num_carbons} out of expected range (18–26) for a C20 precursor"
    
    # Count oxygen atoms.
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    num_oxygens = len(oxygen_atoms)
    if num_oxygens < 3:
        return False, f"Not enough oxygen atoms (found {num_oxygens}; need at least 3)"
    
    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} Da out of expected range (250–500 Da) for nonclassic icosanoids"
    
    # Check for carboxyl functionality.
    # Use three SMARTS patterns: free acid, carboxylate, and ester group.
    carboxyl_smarts = ["C(=O)[OH]", "C(=O)[O-]", "C(=O)O[C]"]
    has_carboxyl = False
    for smarts in carboxyl_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            has_carboxyl = True
            break
    if not has_carboxyl:
        return False, "No carboxyl (acid/carboxylate/ester) group detected"
    
    # Exclude molecules with a 5-membered ring (marker for classic prostanoids).
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            return False, "Contains a 5-membered ring – structure may be a classic prostanoid"
    
    # If all tests pass, we classify as a nonclassic icosanoid.
    return True, "Molecule meets heuristic criteria for nonclassic icosanoid"

# (Optional testing code - remove if not needed)
if __name__ == "__main__":
    test_smiles = [
        "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O",  # (13S)-hydroxy-(14S,15S)-epoxy-(5Z,8Z,11Z)-icosatrienoic acid
        "O[C@@H](CCCC(OC)=O)[C@H](O)/C=C/C=C/C=C\\C=C\\[C@@H](O)CCCCC"   # 5(S),6(R)-Lipoxin A4 methyl ester (should now be detected)
    ]
    for s in test_smiles:
        result, reason = is_nonclassic_icosanoid(s)
        print(result, reason)