"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: A disaccharide, defined as a compound in which two monosaccharides are joined by a glycosidic bond.
A disaccharide is here defined as having two candidate sugar rings (5‐ or 6‐membered rings with exactly one ring oxygen and otherwise carbons)
that also show at least one free hydroxyl substituent and are connected by an oxygen (not in either ring) that bridges a ring carbon from one ring to a ring carbon in the other.
"""

from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is defined as two monosaccharide units (candidate sugar rings) linked via a glycosidic bond.
    
    The algorithm proceeds in a few steps:
      1. Parse the molecule and add explicit hydrogens.
      2. Identify candidate rings. A candidate sugar ring is defined as
         - a 5- or 6-membered ring,
         - containing exactly one oxygen atom (the ring oxygen) and all other members are carbons,
         - and with at least one exocyclic hydroxyl group (an oxygen not in the ring that carries at least one hydrogen).
      3. For a disaccharide we require exactly two candidate rings.
      4. Search for a glycosidic linkage:
         Here the code looks for an oxygen (not belonging to either ring) that is bonded to a carbon in one ring and also to a carbon in the other.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a disaccharide, else False.
        str: A reason explaining the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # add explicit hydrogens to help detect hydroxyl substitution
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each is a tuple of atom indices in a ring

    candidate_rings = []
    for ring in atom_rings:
        # Only consider 5- or 6-membered rings.
        if len(ring) not in (5,6):
            continue
        oxygen_count = 0
        valid_ring = True
        # Check that every atom in the ring is either oxygen or carbon.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            a_num = atom.GetAtomicNum()
            if a_num == 8:
                oxygen_count += 1
            elif a_num == 6:
                # Do not restrict hybridization now.
                pass
            else:
                valid_ring = False
                break
        # Accept only if exactly one oxygen is present.
        if not valid_ring or oxygen_count != 1:
            continue
        
        # Further require that at least one of the ring carbons has a free hydroxyl group.
        # (This checks for an oxygen attached to a ring carbon that is not part of the ring, and that oxygen carries at least one H.)
        hydroxyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # consider only ring carbons
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms in the ring
                if nbr.GetAtomicNum() == 8:
                    # check if this oxygen has a hydrogen neighbor (indicating -OH)
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 1:
                            hydroxyl_found = True
                            break
                    if hydroxyl_found:
                        break
            if hydroxyl_found:
                break
        if not hydroxyl_found:
            continue
        
        candidate_rings.append(set(ring))
        
    # For a disaccharide we expect exactly two candidate sugar rings.
    if len(candidate_rings) != 2:
        return False, f"Expected 2 candidate sugar rings, found {len(candidate_rings)}"
    
    ring1, ring2 = candidate_rings
    
    # Search for a glycosidic linkage.
    # That is, an oxygen atom (that is not part of either candidate ring) bonded to a carbon in ring1 and also to a carbon in ring2.
    glyco_link_found = False
    # Loop over all atoms in ring1 (only carbons) and examine their neighbor oxygens.
    for idx in ring1:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        for nbr in atom.GetNeighbors():
            # Consider only oxygen atoms that are not in any candidate ring.
            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in (ring1.union(ring2)):
                # Check if this oxygen also bonds to a carbon from ring2.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == idx:
                        continue
                    if nbr2.GetIdx() in ring2 and nbr2.GetAtomicNum() == 6:
                        glyco_link_found = True
                        break
                if glyco_link_found:
                    break
        if glyco_link_found:
            break
    # Also check the reverse direction, if not already found.
    if not glyco_link_found:
        for idx in ring2:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in (ring1.union(ring2)):
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == idx:
                            continue
                        if nbr2.GetIdx() in ring1 and nbr2.GetAtomicNum() == 6:
                            glyco_link_found = True
                            break
                    if glyco_link_found:
                        break
            if glyco_link_found:
                break

    if glyco_link_found:
        return True, "Molecule contains two candidate sugar rings with sufficient hydroxyl substitution linked by an exocyclic oxygen, consistent with a disaccharide"
    else:
        return False, "No glycosidic linkage (bridging oxygen connecting ring carbons of the two candidate rings) found"

# Example usage:
# Uncomment one or more of these lines to test with sample disaccharide SMILES.
# ex_smiles = "O([C@@H]1[C@@H](O)[C@@H](O)[C@H](O[C@@H]1O)CO)[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)CO"  # alpha-D-Manp-(1->2)-alpha-D-Galp
# print(is_disaccharide(ex_smiles))