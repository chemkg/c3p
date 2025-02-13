"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophyll
Definition: A subclass of carotenoids consisting of the oxygenated carotenes.
Xanthophylls are characterized by a long conjugated polyene chain (often with at least four consecutive double bonds), an extensive carbon skeleton (typically ≥30 carbon atoms), ring(s) and at least one oxygen-based functional group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    
    A xanthophyll is an oxygenated carotenoid and should have:
    - A long conjugated polyene chain (here, we require a chain containing at least four consecutive C=C bonds).
    - At least one oxygen atom present.
    - A large carbon skeleton (here, at least 30 carbon atoms).
    - At least one ring (since many xanthophylls possess cyclic endgroups).
    - A relatively high molecular weight.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a xanthophyll, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for oxygen atoms (the molecule must be oxygenated)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; not oxygenated"
    
    # Check for a sufficiently large carbon skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Too few carbon atoms ({c_count}); expected at least 30 for a carotenoid backbone"
    
    # Check that the molecule has at least one ring (common in most cyclic carotenoids/xanthophylls)
    if len(mol.GetRingInfo().AtomRings()) < 1:
        return False, "No ring structures found; typically xanthophylls contain at least one cyclic endgroup"
    
    # Look for a conjugated polyene chain.
    # We require at least a chain with 4 consecutive double bonds.
    # The SMARTS below looks for a pattern: C=C-C=C-C=C-C=C
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C-C=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No long conjugated polyene chain (≥4 consecutive C=C bonds) found"
    
    # Check that the molecular weight is high enough (xanthophylls are large molecules)
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 400:
        return False, f"Molecular weight too low ({mw:.1f}); expected at least 400 Da for a xanthophyll"
    
    # All criteria met: classify this molecule as a xanthophyll.
    return True, "Contains a long conjugated polyene chain with ≥4 C=C bonds, sufficient carbon scaffold with ring(s), and oxygen functionalities consistent with xanthophylls"

# Example usage (remove or comment out in production code):
if __name__ == '__main__':
    test_examples = [
        # True positive example (one of the provided ones):
        "O1C2(C1(CCCC2(C)C)C)\\C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C/C=C(\\C3OC4(C(C(CCC4)(C)C)=C3)C)/C)\\C",
        # False positive example (one that was wrongly classified before):
        "OC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)[H])C)[H])C"
    ]
    for sm in test_examples:
        result, reason = is_xanthophyll(sm)
        print("SMILES:", sm)
        print("Result:", result, "|", reason)
        print("-" * 80)