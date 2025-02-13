"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: Hexose
Defined as: Any six‐carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).
This implementation checks if the molecule has exactly 6 carbons and then looks for a carbonyl group or a cyclic sugar ring.
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    
    A hexose is defined as any six‐carbon monosaccharide that in its linear form contains either
    an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).
    This function checks:
      - The molecule contains exactly 6 carbon atoms.
      - If acyclic, it has an aldehyde or ketone moiety.
      - If cyclic, it contains a five‐membered (furanose) or six‐membered (pyranose) ring that includes one oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a hexose, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 6:
        return False, f"Molecule contains {carbon_count} carbon atoms; hexoses should have exactly 6 carbons"
    
    # Define SMARTS patterns for aldehyde and ketone functional groups.
    # Aldehyde: terminal carbonyl with one hydrogen ([CX3H1](=O))
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    # Ketone: carbonyl with carbons on both sides ([#6][CX3](=O)[#6])
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    
    # Check for acyclic carbonyl patterns (open-chain forms)
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Aldehyde group detected; consistent with an aldohexose"
    if mol.HasSubstructMatch(ketone_pattern):
        return True, "Ketone group detected; consistent with a ketohexose"
    
    # If no clear carbonyl is detected then the sugar might be in its cyclic form.
    # Many cyclic hexoses (pyranose or furanose forms) have a ring in which exactly one oxygen is present.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        for ring in ring_info.AtomRings():
            if len(ring) in [5, 6]:
                # Count the number of oxygen atoms in the ring
                oxygens_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                if oxygens_in_ring == 1:
                    return True, f"Cyclic sugar structure detected (ring size {len(ring)} with one oxygen); typical of hexoses"
    
    return False, "Molecule does not match hexose criteria (no appropriate carbonyl found and ring structure not clearly sugar-like)"