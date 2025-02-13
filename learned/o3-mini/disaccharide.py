"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: A disaccharide, defined as a compound in which two monosaccharides are joined by a glycosidic bond.
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is defined as two monosaccharide units (sugar rings) linked via a glycosidic bond.
    
    The algorithm proceeds in two steps:
      1. Identify candidate sugar rings: saturated rings of size 5 or 6 that contain exactly one ring oxygen and 
         otherwise consist of sp3 carbons.
      2. Identify a glycosidic linkage: an oxygen atom (typically not in the rings) that is bonded to one atom 
         in each candidate sugar ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a disaccharide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()    # tuple of tuples; each subtuple is a set of atom indices in a ring

    candidate_rings = []
    # Loop over rings and filter for candidate sugar rings:
    #  - ring size of 5 or 6 
    #  - exactly one oxygen atom within the ring
    #  - all other atoms are carbons in sp3 hybridization
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue
        oxy_count = 0
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                oxy_count += 1
            elif atomic_num == 6:
                # Check that the carbon is sp3-hybridized.
                if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    valid_ring = False
                    break
            else:
                # If any other element appears, skip this ring.
                valid_ring = False
                break
        if valid_ring and oxy_count == 1:
            candidate_rings.append(set(ring))
    
    # For a disaccharide, we expect exactly 2 sugar rings.
    if len(candidate_rings) != 2:
        return False, f"Expected 2 candidate sugar rings, found {len(candidate_rings)}"
    
    ring1, ring2 = candidate_rings

    # Now search for a glycosidic linkage. In many disaccharides the two ring units are linked via an oxygen.
    # Here we look for an oxygen atom (not part of either ring) whose neighbors include at least one atom from ring1 
    # and one atom from ring2.
    glyco_link_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # only interested in oxygen
        # Exclude oxygen if it is part of one of the candidate rings (these are ring oxygens).
        if atom.GetIdx() in ring1 or atom.GetIdx() in ring2:
            continue
        # Get all neighbors of this oxygen.
        neighbor_rings = [False, False]  # [is there neighbor in ring1, ring2]
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in ring1:
                neighbor_rings[0] = True
            if n_idx in ring2:
                neighbor_rings[1] = True
        if all(neighbor_rings):
            glyco_link_found = True
            break

    if not glyco_link_found:
        return False, "No glycosidic bond (bridging oxygen connecting atoms from both sugar rings) found"
    
    return True, "Molecule contains two sugar rings linked by an oxygen glycosidic bond, consistent with a disaccharide"

# Example usage:
# Uncomment the following lines to test the function with a provided SMILES string.
# Example disaccharide: alpha-D-Manp-(1->2)-alpha-D-Galp
# smiles_example = "O([C@@H]1[C@@H](O)[C@@H](O)[C@H](O[C@@H]1O)CO)[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)CO"
# print(is_disaccharide(smiles_example))