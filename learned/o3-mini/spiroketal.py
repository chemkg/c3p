"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: Spiroketal 
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.
Improvement over previous version: 
  - Ensures candidate spiro center is sp³-hybridized and not part of a carbonyl.
  - Only counts oxygen neighbors that belong to a ring together with the candidate.
  - Confirms that for at least one pair of rings the candidate is their only shared atom.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is defined as a cyclic ketal where the ketal carbon (spiro center)
    is the only common atom between two rings. To be a ketal the spiro carbon:
      - must be sp³-hybridized and not be a carbonyl (C=O);
      - must have at least two oxygen neighbors, and those oxygen neighbors
        should be incorporated in at least one ring together with the spiro carbon.
    
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
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Tuple of tuples of atom indices
    
    # Build a dictionary mapping each atom to the list of rings (as sets) that include it.
    atom_to_rings = {}
    for ring in atom_rings:
        ring_set = set(ring)
        for idx in ring:
            atom_to_rings.setdefault(idx, []).append(ring_set)
    
    # Check each carbon atom as candidate spiroketal center.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # must be carbon
        
        idx = atom.GetIdx()
        # Must be in at least two rings
        if idx not in atom_to_rings or len(atom_to_rings[idx]) < 2:
            continue
        
        # Check that the candidate carbon is sp3 hybridized.
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        
        # Exclude if it is already part of a carbonyl: 
        # check if any double bond to oxygen exists.
        has_double_bond_to_O = False
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    has_double_bond_to_O = True
                    break
        if has_double_bond_to_O:
            continue
        
        rings_of_atom = atom_to_rings[idx]
        # Look for at least one pair of rings that share only the candidate carbon.
        spiro_pair_found = False
        for i in range(len(rings_of_atom)):
            for j in range(i+1, len(rings_of_atom)):
                if rings_of_atom[i].intersection(rings_of_atom[j]) == {idx}:
                    spiro_pair_found = True
                    break
            if spiro_pair_found:
                break
        if not spiro_pair_found:
            continue
        
        # For ketal functionality, the candidate carbon should have at least two oxygen neighbors.
        # But we require that these oxygens are incorporated in at least one ring containing the candidate.
        oxy_in_ring_count = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 8:
                continue
            nbr_idx = nbr.GetIdx()
            # Is the oxygen neighbor part of any ring that also contains the candidate?
            if nbr_idx in atom_to_rings:
                # Check in each ring that contains the candidate if this neighbor is also present:
                in_shared_ring = False
                for ring in atom_to_rings[idx]:
                    if nbr_idx in ring:
                        in_shared_ring = True
                        break
                if in_shared_ring:
                    oxy_in_ring_count += 1

        if oxy_in_ring_count < 2:
            return False, f"Found spiro center at atom index {idx} but it only has {oxy_in_ring_count} oxygen neighbors in rings."
        
        # Passed all tests: candidate spiro center shows the traits of a spiroketal.
        return True, f"Found spiroketal center at atom index {idx} shared by two rings with at least 2 ring-incorporated oxygen substituents."
    
    return False, "No spiroketal center identified (no carbon found that is the sole common atom between two rings with proper ketal oxygen substituents)."

# Uncomment below to test examples:
# test_smiles = "O=C1O[C@@H]2C[C@@]3(O[C@H](C(=CC)C)[C@@H](C)C[C@H]3O)O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C(C(=O)C4)C)O)CO)C"  # VM-44864
# print(is_spiroketal(test_smiles))