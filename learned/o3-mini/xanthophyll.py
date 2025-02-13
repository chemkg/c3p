"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophyll
Definition: A subclass of carotenoids consisting of the oxygenated carotenes.
Xanthophylls are characterized by a long conjugated polyene chain, an extensive carbon skeleton,
and oxygen-based functional groups. Many have cyclic endgroups, but some acyclic xanthophylls exist.
For classification, we require:
  - At least one oxygen atom.
  - A sufficiently large carbon backbone (≥30 carbons).
  - A long conjugated polyene chain: if a ring is present we require at least 4 consecutive C=C bonds;
    if no ring is present (acyclic) we require a longer conjugated chain (≥5 consecutive C=C bonds).
  - A relatively high molecular weight (≥400 Da).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    
    A xanthophyll is defined as an oxygenated carotenoid. Here, we enforce that the molecule:
      - Has at least one oxygen atom.
      - Contains at least 30 carbon atoms.
      - Possesses a long conjugated polyene chain.
         * If at least one ring exists, we require a pattern of 4 consecutive C=C bonds.
         * If no ring exists (acyclic candidate), we require a pattern of 5 consecutive C=C bonds.
      - Has a molecular weight of at least 400 Da.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a xanthophyll, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check oxygenation
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; not an oxygenated carotene"
    
    # Check carbon scaffold; require at least 30 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Too few carbon atoms ({c_count}); expected at least 30 for a carotenoid backbone"
    
    # Molecular weight check (xanthophylls are large molecules)
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 400:
        return False, f"Molecular weight too low ({mw:.1f} Da); expected ≥400 Da for a xanthophyll"
    
    # Check for rings.
    rings = mol.GetRingInfo().AtomRings()
    has_ring = len(rings) > 0
    
    # Depending on whether a ring exists, choose a SMARTS pattern for a conjugated polyene chain.
    # (Note: these SMARTS patterns are simplified; more elaborate conjugation detection is possible.)
    if has_ring:
        # For cyclic xanthophylls, require at least 4 consecutive C=C bonds.
        polyene_smarts = "C=C-C=C-C=C-C=C"
        chain_info = "a conjugated chain with at least 4 consecutive C=C bonds"
    else:
        # For acyclic xanthophylls, require a longer chain (here: 5 consecutive C=C bonds).
        polyene_smarts = "C=C-C=C-C=C-C=C-C=C"
        chain_info = "a conjugated chain with at least 5 consecutive C=C bonds (typical for acyclic carotenoids)"
    
    polyene_pattern = Chem.MolFromSmarts(polyene_smarts)
    if polyene_pattern is None:
        return False, "Internal error: could not generate polyene SMARTS pattern"
    
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, f"No long conjugated polyene chain found; expected {chain_info}"
    
    # All conditions met.
    ring_info = "with at least one ring" if has_ring else "without ring structures (acyclic)"
    return True, (f"Contains {ring_info}, a long conjugated polyene chain, sufficient carbon scaffold "
                  f"({c_count} C atoms) and oxygen functionalities (O count: {o_count}), "
                  f"with molecular weight {mw:.1f} Da consistent with xanthophylls")

# Example usage (for testing); remove or comment out before production deployment:
if __name__ == '__main__':
    test_examples = [
        # A cyclic example: (5R,5'R,6S,8'R)-Luteochrome
        "O1C2(C1(CCCC2(C)C)\\C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C/C=C(\\C3OC4(C(C(CCC4)(C)C)=C3)C)/C)\\C",
        # An acyclic xanthophyll candidate: spheroidene
        "COC(C)(C)C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)CC\\C=C(/C)CCC=C(C)C",
    ]
    for sm in test_examples:
        result, reason = is_xanthophyll(sm)
        print("SMILES:", sm)
        print("Result:", result, "|", reason)
        print("-" * 80)