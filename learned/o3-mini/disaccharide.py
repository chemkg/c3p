"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: A disaccharide, defined as a compound in which two monosaccharides are joined by a glycosidic bond.
A disaccharide is here defined as having two candidate sugar rings (5‐ or 6‐membered rings with one ring oxygen
and otherwise sp3 carbons) that are connected by a bridging oxygen.
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is defined as two monosaccharide units (candidate sugar rings) linked via a glycosidic bond.
    
    The algorithm proceeds in two main steps:
      1. Identify candidate sugar rings. Here a candidate ring is a 5‐ or 6‐membered ring that contains exactly
         one oxygen atom (as typically found in pyranose or furanose rings) and all other members are carbons that are sp3‐hybridized.
      2. Identify a glycosidic linkage. Instead of only looking at oxygens “outside” the rings we scan every candidate ring carbon’s neighbors.
         If an oxygen attached to a ring carbon is also bonded to a carbon that belongs to the other candidate ring, we count that as a glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a disaccharide, else False.
        str: A reason explaining the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples (each tuple is a set of atom indices in a ring)
    
    candidate_rings = []
    # Loop over rings to filter candidate sugar rings:
    #   - only proceed if ring size is 5 or 6.
    #   - require exactly one oxygen atom within the ring.
    #   - require that all the non‐oxygen atoms are carbons and sp3 hybridized.
    for ring in atom_rings:
        if len(ring) not in (5,6):
            continue
        oxygen_count = 0
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:  # oxygen in the ring
                oxygen_count += 1
            elif atomic_num == 6:  # carbon atom
                # Check that it is sp3 hybridized.
                if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    valid_ring = False
                    break
            else:
                # If any other element appears in the ring, reject it.
                valid_ring = False
                break
        if valid_ring and oxygen_count == 1:
            candidate_rings.append(set(ring))
    
    # For a disaccharide we expect exactly 2 candidate sugar rings.
    if len(candidate_rings) != 2:
        return False, f"Expected 2 candidate sugar rings, found {len(candidate_rings)}"
    
    ring1, ring2 = candidate_rings

    # Now search for a glycosidic linkage.
    # Instead of scanning over all (non‐ring) oxygens, we look at each carbon in ring1 and ring2.
    # For each ring carbon, we check if any of its oxygen neighbors is also bonded to a carbon in the opposite ring.
    glyco_link_found = False

    # Check from ring1 to ring2.
    for idx in ring1:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider carbons (typically the anomeric carbon is in the ring).
        if atom.GetAtomicNum() != 6:
            continue
        for neighbor in atom.GetNeighbors():
            # Look for an oxygen.
            if neighbor.GetAtomicNum() == 8:
                # Examine all neighbors of this oxygen (except our original carbon).
                for neigh2 in neighbor.GetNeighbors():
                    if neigh2.GetIdx() == idx:
                        continue
                    if neigh2.GetIdx() in ring2:
                        glyco_link_found = True
                        break
                if glyco_link_found:
                    break
        if glyco_link_found:
            break
    
    # Also check the reverse direction:
    if not glyco_link_found:
        for idx in ring2:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    for neigh2 in neighbor.GetNeighbors():
                        if neigh2.GetIdx() == idx:
                            continue
                        if neigh2.GetIdx() in ring1:
                            glyco_link_found = True
                            break
                    if glyco_link_found:
                        break
            if glyco_link_found:
                break

    if not glyco_link_found:
        return False, "No glycosidic bond (bridging oxygen linking atoms from both candidate sugar rings) found"
    else:
        return True, "Molecule contains two sugar rings linked by a glycosidic oxygen, consistent with a disaccharide"

# Example usage:
# Uncomment the following lines to test the function with one SMILES.
# Example disaccharide: alpha-D-Manp-(1->2)-alpha-D-Galp
# smiles_example = "O([C@@H]1[C@@H](O)[C@@H](O)[C@H](O[C@@H]1O)CO)[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)CO"
# print(is_disaccharide(smiles_example))