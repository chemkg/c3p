"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: nonclassic icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
other than the classic icosanoids (the leukotrienes and the prostanoids).
Heuristic rules used:
  1. Molecule must be a valid structure.
  2. It should have at least ~20 carbon atoms (allowing some wiggle room for additional substituents).
  3. It contains at least one carboxylic acid group.
  4. It must be oxygenated (contains several oxygen atoms).
  5. Molecules containing a 5‐membered ring are assumed to be classic prostanoids and are excluded.
Note: This is a heuristic approach and may not cover all edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    
    Heuristic criteria:
      - Valid molecule.
      - Total number of carbon atoms is at least 20.
      - Contains at least one carboxylic acid group (–C(=O)O).
      - Has a significant number of oxygen atoms (>=3).
      - Does not contain a 5-membered ring, a marker for prostanoids.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a nonclassic icosanoid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbon atoms (atomic number 6)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)
    if num_carbons < 20:
        return False, f"Not enough carbons (found {num_carbons}, need at least 20)"
    
    # Count total oxygen atoms (atomic number 8)
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    num_oxygens = len(oxygen_atoms)
    if num_oxygens < 3:
        return False, f"Not enough oxygen atoms (found {num_oxygens}, need at least 3)"
    
    # Check for carboxylic acid group using a SMARTS: carboxyl acid = C(=O)[OH] or its deprotonated form.
    # Here we use two SMARTS to capture acid and carboxylate.
    acid_smarts_list = ["C(=O)[OH]", "C(=O)[O-]"]
    acid_found = False
    for sm in acid_smarts_list:
        acid_pattern = Chem.MolFromSmarts(sm)
        if mol.HasSubstructMatch(acid_pattern):
            acid_found = True
            break
    if not acid_found:
        return False, "No carboxylic acid group detected"
    
    # Check for classic prostanoid feature: a cyclopentane ring.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    for ring in atom_rings:
        if len(ring) == 5:
            return False, "Contains a 5-membered ring - classic prostanoid structure suspected"
            
    # Additional heuristics could be implemented here to exclude leukotrienes,
    # but given the complexity we assume that molecules passing the above criteria
    # are likely nonclassic icosanoids.
    
    return True, "Molecule meets heuristic criteria for nonclassic icosanoid"