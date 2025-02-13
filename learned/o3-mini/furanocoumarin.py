"""
Classifies: CHEBI:24128 furanocoumarin
"""
#!/usr/bin/env python
"""
Classifies: Any furochromene that consists of a furan ring fused with a coumarin.
The function is_furanocoumarin uses SMARTS patterns for the coumarin (benzopyrone) core
and then searches the ring information for aromatic five‐membered rings with one oxygen atom.
Fusion is defined as the requirement that at least one candidate furan ring shares two or more atoms with one coumarin match.
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a coumarin (benzopyrone) having a fused furan ring.
    Fusion is confirmed if a candidate coumarin match and a candidate furan ring share at least 2 atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a furanocoumarin, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the coumarin (benzopyrone) core.
    coumarin_smarts_list = [
        "c1ccc2oc(=O)cc2c1",   # common representation
        "c1ccc2c(c1)oc(=O)cc2"  # a slight variation
    ]
    
    coumarin_matches = []
    for smarts in coumarin_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            # Add each coumarin substructure match
            coumarin_matches.extend(matches)
    
    if not coumarin_matches:
        return False, "Coumarin (benzopyrone) core not found"
    
    # Use ring information to identify candidate furan rings.
    # Candidate furan ring: 5 atoms, exactly one oxygen, the atoms must be aromatic.
    ring_info = mol.GetRingInfo()
    candidate_furan_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        # Count oxygens and check aromaticity.
        o_count = 0
        aromatic_flag = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                o_count += 1
            if not atom.GetIsAromatic():
                aromatic_flag = False
                break
        # Candidate if exactly one oxygen, four carbons, and all atoms aromatic.
        if o_count == 1 and aromatic_flag:
            # Additionally check that the rest are carbons.
            carbon_flag = all(mol.GetAtomWithIdx(idx).GetAtomicNum() in (6, 8) for idx in ring)
            if carbon_flag:
                candidate_furan_rings.append(set(ring))
    
    if not candidate_furan_rings:
        return False, "Aromatic furan ring (5-membered ring with one oxygen) not found"
    
    # Check for fusion: require at least one coumarin match and one candidate furan ring share two or more atoms.
    for c_match in coumarin_matches:
        set_c = set(c_match)
        for f_ring in candidate_furan_rings:
            common_atoms = set_c.intersection(f_ring)
            if len(common_atoms) >= 2:
                return True, "Contains a coumarin (benzopyrone) core fused with an aromatic furan ring (≥2 shared atoms)"
    
    return False, "Coumarin core and furan ring not properly fused (shared atoms < 2)"

# Example usage:
if __name__ == "__main__":
    # Test with a known furanocoumarin (Isobergapten)
    test_smiles = "COc1cc2occc2c2oc(=O)ccc12"
    result, reason = is_furanocoumarin(test_smiles)
    print(result, ":", reason)