"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: Anilide 
Definition: Anilide is any aromatic amide obtained by acylation of aniline.
That is, the amide (C(=O)-N) group must have the acyl carbon attached to nitrogen and
the nitrogen must be directly bonded (as a substituent) to an aromatic (benzene-type) ring.
Additionally, we require that aside from the carbonyl carbon, the amide nitrogen (which
originally belonged to aniline) has no extra heavy substituents.
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    Anilides are aromatic amides formed by acylation of aniline,
    meaning the amide nitrogen (N in C(=O)N) must be directly bound to a benzene ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is an anilide, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Define an amide bond SMARTS pattern: a carbonyl carbon attached to a nitrogen.
    # This is a simple match that will match any C(=O)N fragment.
    amide_pattern = Chem.MolFromSmarts("[#6](=O)[#7]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found"
    
    # Iterate through each identified amide bond match
    # The match order: (carbonyl carbon, oxygen, nitrogen)
    for match in amide_matches:
        carbonyl_c_idx = match[0]
        nitrogen_idx = match[2]
        n_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Get all heavy (non-hydrogen) neighbors of the amide nitrogen.
        heavy_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        
        # For anilides derived from aniline, after acylation the N should have exactly two heavy neighbors:
        # one is the carbonyl carbon and the other is the aromatic (benzene) carbon.
        if len(heavy_neighbors) != 2:
            # Not matching the expected substitution pattern for aniline (could be extra substituents)
            continue
        
        # Identify the neighbor that is not the carbonyl carbon.
        non_carbonyl_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetIdx() != carbonyl_c_idx]
        if not non_carbonyl_neighbors:
            continue
        
        # There should be exactly one non-carbonyl heavy neighbor.
        candidate = non_carbonyl_neighbors[0]
        # Check that the candidate is an aromatic carbon.
        if candidate.GetAtomicNum() != 6 or not candidate.GetIsAromatic():
            continue
        
        # Further require that this aromatic neighbor is part of a six‐membered ring.
        in_six_membered_ring = False
        for ring in ring_info:
            if candidate.GetIdx() in ring and len(ring) == 6:
                in_six_membered_ring = True
                break
        if not in_six_membered_ring:
            continue
        
        # If we reach here, then the amide bond appears to be from aniline acylation.
        return True, "Found an amide bond with the N attached to a benzene ring (6-membered aromatic) indicative of anilide structure."
    
    return False, "No amide bond with the nitrogen bound to a benzene ring (as expected in anilides) was found."
    
# Uncomment below for simple testing:
# test_smiles = "CCCOCCN(C(=O)CCl)c1c(CC)cccc1CC"  # Pretilachlor from examples
# print(is_anilide(test_smiles))