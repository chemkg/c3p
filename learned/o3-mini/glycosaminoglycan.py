"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: Glycosaminoglycan
Definition: Any polysaccharide containing a substantial proportion of aminomonosaccharide residues.
Heuristic: Identify candidate sugar rings (5 or 6‐membered rings containing one oxygen)
           then check if the ring has a substituent that is an amino group.
           The molecule is classified as a glycosaminoglycan if it has at least 2 candidate sugar rings
           and at least 30% of these show amino substitution.
"""

from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    
    Glycosaminoglycans are defined as polysaccharides (multiple repeating carbohydrate units)
    that contain a substantial proportion of aminomonosaccharide residues. As a heuristic,
    we look for candidate sugar rings -- defined here as 5- or 6-membered rings with exactly one ring oxygen --
    and among these, we consider a ring an aminomonosaccharide if at least one of its ring atoms
    bears a substituent nitrogen (i.e. not part of the ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as glycosaminoglycan, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_sugar_rings = []
    aminomonosaccharide_count = 0
    
    # Loop over all rings in the molecule
    for ring in atom_rings:
        # Consider rings of size 5 or 6 as candidate sugar rings.
        if len(ring) not in (5, 6):
            continue
        
        # Count how many atoms in the ring are oxygen.
        oxygen_count = 0
        other_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            else:
                other_count += 1
        
        # For many sugar rings (furanose or pyranose), we expect exactly one ring oxygen.
        if oxygen_count != 1:
            continue
        
        # Passed the simple ring criteria; add to candidate list.
        candidate_sugar_rings.append(ring)
        
        # Check for an amino substituent attached to any ring atom.
        found_amino = False
        ring_set = set(ring)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look for substituents that are outside the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                # If a nitrogen is found attached, consider this ring as an aminomonosaccharide.
                if nbr.GetAtomicNum() == 7:
                    found_amino = True
                    break
            if found_amino:
                break
        
        if found_amino:
            aminomonosaccharide_count += 1

    # For a polysaccharide we expect at least 2 sugar rings.
    total_sugars = len(candidate_sugar_rings)
    if total_sugars < 2:
        return False, f"Only {total_sugars} candidate sugar ring(s) found; need a polysaccharide chain."

    # Now check if a substantial proportion of the sugar rings are aminomonosaccharides
    fraction_amino = aminomonosaccharide_count / total_sugars
    if fraction_amino < 0.3:
        return False, (f"Only {aminomonosaccharide_count} out of {total_sugars} sugar rings "
                       f"({fraction_amino*100:.1f}%) have an amino substituent; "
                       "not enough aminomonosaccharide residues.")

    # If criteria are met, classify as glycosaminoglycan.
    return True, (f"Detected {total_sugars} candidate sugar rings, with {aminomonosaccharide_count} "
                  f"aminosugar(s) ({fraction_amino*100:.1f}%), consistent with a glycosaminoglycan.")

# Example usage:
if __name__ == "__main__":
    test_smiles = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2OC(CO)C(O)C(O)C2O)[C@@H](O)[C@@H]1O"  # A simplified sugar example
    result, reason = is_glycosaminoglycan(test_smiles)
    print(result, reason)