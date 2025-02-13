"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: Spiroketal
Definition: A cyclic ketal in which the ketal (spiro) carbon is the only common atom of two rings.
Improvements:
 - The candidate spiro center (carbon) is required to be sp3‐hybridized and not part of a carbonyl.
 - The candidate carbon must have exactly two oxygen substituents via single bonds.
 - Instead of only checking that there is a pair of rings that share the candidate (the minimum definition
   of a spiro center), we ensure that for at least one pair of rings that only share the candidate,
   one oxygen substituent is present exclusively in one ring and the other oxygen exclusively in the other.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    
    A spiroketal is defined as a cyclic ketal in which:
     - The ketal (spiro) carbon is sp3-hybridized and not part of a carbonyl.
     - The spiro carbon has exactly two oxygen (OR) substituents via single bonds.
     - There exists at least one pair of rings that share only the candidate carbon, and in these rings,
       one oxygen substituent is present in one ring and the other in the other ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a spiroketal, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information (each ring is a tuple of atom indices).
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Build a mapping: for each atom index, list of rings (as sets) that include that atom.
    atom_to_rings = {}
    for ring in atom_rings:
        ring_set = set(ring)
        for idx in ring:
            atom_to_rings.setdefault(idx, []).append(ring_set)
    
    # Loop over atoms looking for candidate spiro centers.
    for atom in mol.GetAtoms():
        # Only consider carbon atoms as candidate spiro centers.
        if atom.GetAtomicNum() != 6:
            continue
        idx = atom.GetIdx()
        
        # Candidate must be in at least two rings.
        if idx not in atom_to_rings or len(atom_to_rings[idx]) < 2:
            continue
        
        # Candidate must be sp3-hybridized.
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        
        # Exclude candidate if it is part of any carbonyl (double bond to oxygen).
        is_carbonyl = False
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    is_carbonyl = True
                    break
        if is_carbonyl:
            continue
        
        # For ketal functionality: the candidate must have exactly two oxygen neighbors 
        # connected by single bonds.
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
            if bond.GetBondType() != rdchem.BondType.SINGLE:
                continue
            if nbr.GetAtomicNum() == 8:
                oxy_neighbors.append(nbr)
        if len(oxy_neighbors) != 2:
            continue
        
        # Get the list of rings that candidate belongs to.
        rings_of_candidate = atom_to_rings[idx]
        
        # Now, check for a pair of rings that share ONLY the candidate and where we can assign
        # one distinct oxygen neighbor to each ring.
        spiro_valid = False
        for i in range(len(rings_of_candidate)):
            for j in range(i+1, len(rings_of_candidate)):
                ring_i = rings_of_candidate[i]
                ring_j = rings_of_candidate[j]
                # Check that the two rings share only the candidate atom.
                if ring_i.intersection(ring_j) != {idx}:
                    continue
                
                # For each ring, count how many of the candidate's oxygen neighbors are present.
                oxy_in_ring_i = [oxy for oxy in oxy_neighbors if oxy.GetIdx() in ring_i]
                oxy_in_ring_j = [oxy for oxy in oxy_neighbors if oxy.GetIdx() in ring_j]
                
                # In a proper spiroketal, one ring should contain exactly one oxygen neighbor and the other the other.
                if len(oxy_in_ring_i) == 1 and len(oxy_in_ring_j) == 1:
                    # Extra check: ensure they are distinct oxygen atoms.
                    if oxy_in_ring_i[0].GetIdx() != oxy_in_ring_j[0].GetIdx():
                        spiro_valid = True
                        break
            if spiro_valid:
                break
        
        if spiro_valid:
            return True, (f"Found spiroketal center at atom index {idx} "
                          "with one oxygen substituent in one ring and the other in a different ring "
                          "that share only the candidate carbon.")
    
    return False, "No spiroketal center identified that meets the defined ketal and spiro topology criteria."

# For testing (uncomment if needed):
# test_smiles = "O=C1O[C@@H]2C[C@@]3(O[C@H](C(=CC)C)[C@@H](C)C[C@H]3O)O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C(C(=O)C4)C)O)CO)C" 
# print(is_spiroketal(test_smiles))