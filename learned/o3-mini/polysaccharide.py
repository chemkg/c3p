"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
It is typically used for structures with more than ten monosaccharide residues.
The detection here is based on the identification of sugar rings (furanose or pyranose rings)
defined as rings of size 5 (1 oxygen, 4 carbons) or size 6 (1 oxygen, 5 carbons).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    
    Criteria:
      - The molecule must be valid.
      - It must contain >10 monosaccharide residues. 
        For simplicity, this function identifies monosaccharide residues as rings of size 5 or 6 that
        contain exactly one oxygen and the remaining atoms as carbons. (Furanose: 5-membered, Pyranose: 6-membered)
      - A molecular weight check is also applied (e.g., >1000 Da) to confirm the structure is biomacromolecular.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a polysaccharide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information
    ring_info = mol.GetRingInfo()
    sugar_ring_count = 0
    
    # Loop over all rings in the molecule
    for ring in ring_info.AtomRings():
        # Only consider 5- or 6-membered rings
        if len(ring) in (5, 6):
            # Count oxygen atoms and carbon atoms in the ring
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            carbon_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # For a typical sugar unit:
            # 5-membered ring: 1 oxygen and 4 carbons.
            # 6-membered ring: 1 oxygen and 5 carbons.
            if oxygen_count == 1 and ((len(ring) == 5 and carbon_count == 4) or (len(ring) == 6 and carbon_count == 5)):
                sugar_ring_count += 1
    
    # Polysaccharide requirement: more than 10 sugar rings (i.e., at least 11 monosaccharide residues)
    if sugar_ring_count <= 10:
        return False, f"Found only {sugar_ring_count} sugar ring(s); a polysaccharide requires more than 10 residues."
    
    # Check if the molecular weight is in the expected range (e.g., greater than 1000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a typical polysaccharide."
    
    return True, f"Detected {sugar_ring_count} sugar ring(s) with a molecular weight of {mol_wt:.1f} Da indicative of a polysaccharide."