"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: Spiroketal
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.
Improvements over previous version:
  - Candidate spiro centers are required to be sp³‐hybridized carbon atoms that are not part of a carbonyl.
  - Candidate carbon must have exactly two oxygen substituents (via single bonds) – as expected for a ketal.
  - Only if both oxygen neighbors occur in at least one candidate ring and if there is at least one pair of rings 
    that share exclusively the candidate, the structure is classified as a spiroketal.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    
    A spiroketal is defined as a cyclic ketal where:
      - The ketal (spiro) carbon must be sp3-hybridized and not be part of a carbonyl.
      - The ketal carbon should form exactly two single bonds to oxygen (that is, two OR groups).
      - The spiro center must lie in at least two rings, and there must be at least one pair of rings 
        that intersect (share atoms) only at the candidate carbon.
      - Each oxygen substituent must be part of at least one ring together with the candidate.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a spiroketal, False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Tuple of tuples, each tuple is a ring (atom indices)
    
    # Build a dictionary mapping each atom index to a list of rings (represented as sets) that include that atom.
    atom_to_rings = {}
    for ring in atom_rings:
        ring_set = set(ring)
        for idx in ring:
            atom_to_rings.setdefault(idx, []).append(ring_set)
    
    # Loop over atoms looking for candidate spiroketal centers.
    for atom in mol.GetAtoms():
        # Only consider carbon atoms as the potential spiro center.
        if atom.GetAtomicNum() != 6:
            continue
        
        idx = atom.GetIdx()
        # Candidate must be in at least two rings.
        if idx not in atom_to_rings or len(atom_to_rings[idx]) < 2:
            continue
        
        # Candidate must be sp3 hybridized.
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        
        # Exclude candidate if it is part of any carbonyl (i.e. has any double bond to oxygen).
        is_carbonyl = False
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    is_carbonyl = True
                    break
        if is_carbonyl:
            continue
        
        # For ketal functionality, the candidate should have exactly two oxygen neighbors via single bonds.
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Only consider single bonds.
            if bond.GetBondType() != rdchem.BondType.SINGLE:
                continue
            if nbr.GetAtomicNum() == 8:
                oxy_neighbors.append(nbr)
        if len(oxy_neighbors) != 2:
            # Not exactly two oxygen substituents – skip this candidate.
            continue
        
        # Verify that there exists at least one pair of rings containing the candidate
        # that share no atom other than the candidate.
        rings_of_candidate = atom_to_rings[idx]
        spiro_pair_found = False
        for i in range(len(rings_of_candidate)):
            for j in range(i+1, len(rings_of_candidate)):
                if rings_of_candidate[i].intersection(rings_of_candidate[j]) == {idx}:
                    spiro_pair_found = True
                    break
            if spiro_pair_found:
                break
        if not spiro_pair_found:
            continue
        
        # For each oxygen neighbor, confirm that it is part of at least one ring that also contains the candidate.
        valid_oxy_count = 0
        for oxy in oxy_neighbors:
            oxy_idx = oxy.GetIdx()
            in_shared_ring = False
            for ring in rings_of_candidate:
                if oxy_idx in ring:
                    in_shared_ring = True
                    break
            if in_shared_ring:
                valid_oxy_count += 1
        if valid_oxy_count != 2:
            return False, f"Found candidate spiro center at atom index {idx} but only {valid_oxy_count} oxygen substituents are in rings with it."
        
        # Candidate passed all tests. Report success.
        return True, f"Found spiroketal center at atom index {idx} that is the sole shared atom between two rings and has exactly 2 ring-incorporated oxygen substituents."

    return False, "No spiroketal center identified with the required ketal functionality and spiro ring topology."


# Uncomment the following lines to test examples:
# test_smiles = "O=C1O[C@@H]2C[C@@]3(O[C@H](C(=CC)C)[C@@H](C)C[C@H]3O)O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C(C(=O)C4)C)O)CO)C"  # Example: VM-44864
# print(is_spiroketal(test_smiles))