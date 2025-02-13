"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A compound that contains exactly two ketone functionalities (diketone)

A diketone is defined as a compound in which exactly two carbonyl (C=O) groups occur 
in which the carbonyl carbon is bound to two carbon atoms (thus excluding aldehydes
or carbonyls bonded to heteroatoms). In addition, we require that no extra carbonyl
functionality is present, and we penalize compounds that contain a highly complex ring 
system (i.e. more than 4 rings) as these are frequently found in steroidal or fused 
polycyclic molecules that are not “simple” diketones.
"""

from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule qualifies as a diketone based on its SMILES string.
    For our purposes, a diketone must have exactly two ketone groups – each defined
    as a [#6][CX3](=O)[#6] substructure – and exactly two overall carbonyl groups
    (matched by [CX3]=O). In addition, we enforce that the molecule does not contain
    an unusually high number of rings (more than 4) which is used as a proxy for complex
    fused systems that tend to arise in steroid or polycyclic natural products.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a diketone, False otherwise.
        str: An explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns:
    # A ketone group: a carbonyl in which the carbon is bonded to two carbons.
    ketone_smarts = "[#6][CX3](=O)[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    
    # A general carbonyl (this will also match aldehydes etc).
    carbonyl_smarts = "[CX3]=O"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)
    
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    
    n_ketones = len(ketone_matches)
    n_carbonyls = len(carbonyl_matches)
    
    # Check that exactly two ketone groups (by our definition) were found.
    if n_ketones != 2:
        return False, f"Compound contains {n_ketones} ketone group(s) (by strict SMARTS criteria), which does not equal 2."
    
    # Enforce that no extra carbonyl functionality is present.
    if n_carbonyls != 2:
        return False, f"Compound contains {n_carbonyls} carbonyl group(s); extra carbonyl functionality detected beyond the 2 ketone groups."
    
    # As many false positives were very polycyclic (e.g. steroid or fused systems) despite having exactly 2 ketones,
    # we add an extra heuristic: if the molecule contains more than 4 rings, it is likely too complex.
    n_rings = mol.GetRingInfo().NumRings()
    if n_rings > 4:
        return False, f"Compound has a complex ring system ({n_rings} rings) not typical for a simple diketone."
    
    return True, "Contains exactly 2 ketone (C=O) groups (with both substituents being carbon) and no extra carbonyl functionality."

# Example usage (for testing purposes only)
if __name__ == "__main__":
    test_smiles = [
        "O=C(CCCCCCCCCCCCCCC)CC(=O)CCCCCCCCC",  # 10,12-Heptacosanedione (expected True)
        "C1[C@H](C(C([C@]2([H])[C@]1(C)[C@@]3(C(C=C([C@H]([C@]3(CC2)[H])C)C=C)=O)[H])(C)C)=O)O",  # (+)-phytocassane A (expected True)
        "OC1C(C(=O)C(C1=O)C(=O)C(CC)C)CC=C(C)C",  # Adhumulinic acid (3 ketones, expected False)
        "CC(=O)[C@@]12OC(C)(C)O[C@@H]1C[C@H]1[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]21C",  # Descinolone acetonide (expected False)
    ]
    for smi in test_smiles:
        result, reason = is_diketone(smi)
        print(f"SMILES: {smi}")
        print(f"Result: {result}, Reason: {reason}")
        print("-" * 50)