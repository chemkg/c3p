"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Definition: Any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
A candidate is identified by finding an isolated (non‐fused) benzene ring (six aromatic carbons) and then checking that
exactly 3 of the ring carbons bear a terminal chlorine substituent. In addition, any other substituent attached to the ring
(which is not an aromatic fragment) must be “small” (here, defined as having 2 or fewer heavy atoms).
"""

from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule qualifies as a trichlorobenzene based on its SMILES string.
    
    The algorithm:
      1. Parse and sanitize the molecule.
      2. Identify rings that are candidate benzene fragments:
         - a ring of exactly 6 atoms,
         - all atoms are aromatic carbons.
         - the ring is isolated (not fused with any other ring; i.e. it does not share 2 or more atoms with any other ring).
      3. For each candidate ring, examine each ring atom’s neighbors that are not in the ring:
         - If the neighbor is a chlorine (atomic number 17) and is terminal (degree==1) then count it.
         - Otherwise, perform a DFS starting from that neighbor – but only “count” heavy atoms if no aromatic atom is encountered.
           (If an aromatic atom is reached, we assume this bond is simply the linkage to an aromatic fragment and allow it.)
         - If any substituent (other than Cl) is found to span more than 2 heavy atoms [without quickly “hitting” an aromatic atom],
           then disqualify that candidate benzene ring.
      4. If a candidate ring has exactly 3 valid terminal Cl substituents (and no disqualifying larger non‐Cl substituents),
         return True plus an explanation.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if a candidate ring meeting the criteria is found, False otherwise.
        str: Explanation of the classification decision.
    """
    # Helper DFS that traverses a substituent fragment (avoiding the candidate ring atom that initiated the bond)
    # The DFS returns a tuple (count, aromatic_found) where count is the total number of heavy atoms (atomic number > 1)
    # encountered in the fragment (if no aromatic atom is encountered) and aromatic_found is True if an aromatic atom is reached.
    def dfs_substituent(atom, origin_idx, visited):
        # If we hit an aromatic atom, then we mark that branch as aromatic.
        if atom.GetIsAromatic():
            return (0, True)
        # Count this heavy atom (ignore hydrogens by atomic number <= 1)
        count = 1 if atom.GetAtomicNum() > 1 else 0
        aromatic_found = False
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == origin_idx or nbr.GetIdx() in visited:
                continue
            sub_count, sub_arom = dfs_substituent(nbr, origin_idx, visited)
            if sub_arom:
                aromatic_found = True
            else:
                count += sub_count
            # If already aromatic found in one branch, we can stop further counting as this branch is allowed.
            if aromatic_found:
                # No need to accumulate further heavy atoms.
                break
        return (count, aromatic_found)
    
    # Parse the input SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Loop over all rings in the molecule to find a candidate benzene ring.
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        
        # Check that every atom in the ring is a carbon and is aromatic.
        is_candidate = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_candidate = False
                break
        if not is_candidate:
            continue
            
        # Exclude fused rings: if the candidate ring shares 2 or more atoms with any other ring in the molecule.
        fused = False
        for other_ring in atom_rings:
            if other_ring == ring:
                continue
            if len(set(ring).intersection(other_ring)) >= 2:
                fused = True
                break
        if fused:
            continue
        
        # For the candidate benzene ring, check the substituents.
        cl_count = 0  # count of terminal chlorine substituents on the ring.
        disqualified = False
        
        # Iterate over each atom in the candidate ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Skip neighbors that are part of the ring.
                if nbr.GetIdx() in ring:
                    continue
                # If the neighbor is an aromatic atom, then we assume it is an attached aromatic fragment.
                if nbr.GetIsAromatic():
                    continue
                # If the neighbor is chlorine (atomic #17) and is terminal (degree == 1), count it.
                if nbr.GetAtomicNum() == 17 and nbr.GetDegree() == 1:
                    cl_count += 1
                else:
                    # For other substituents, assess its "size" using DFS.
                    # Create a new visited set for each substituent branch.
                    visited = set()
                    count, arom_found = dfs_substituent(nbr, atom.GetIdx(), visited)
                    # If an aromatic fragment is encountered in this substituent branch,
                    # then we treat it as not disqualifying (assumed to be a simple aromatic extension).
                    if not arom_found and count > 2:
                        disqualified = True
                        # We can break out if one substituent is too extensive.
                        break
            if disqualified:
                break
        
        # We accept this candidate ring only if there are exactly 3 terminal Cl substituents and no disqualifying substituents.
        if not disqualified and cl_count == 3:
            return True, ("Found an isolated benzene ring with exactly three terminal chlorine substituents "
                          "and only small (or aromatic) additional substituents")
    
    return False, ("No isolated benzene ring with exactly three terminal chlorine substituents "
                   "(with acceptable small or aromatic substituents) found")


# (Optional) Example testing
if __name__ == '__main__':
    test_smiles = [
        # Expected True examples:
        "Oc1c(Cl)cc(c(Cl)c1Cl)-c1ccc(Cl)c(Cl)c1Cl",  # 2,2',3,3',4',5-Hexachloro-4-biphenylol
        "Clc1cccc(Cl)c1Cl",  # 1,2,3-trichlorobenzene
        # Some examples that should be False:
        "C=1(C(=C(C=C(C1Cl)O)Cl)[O-])Cl",  # 2,3,6-trichloro-4-hydroxyphenolate (false positive in previous attempt)
        "CON([C@@H](C)Cc1c(Cl)cc(Cl)cc1Cl)C(=O)c1cn(C)nc1C(F)F"  # (S)-pydiflumetofen (false negative in previous attempt)
    ]
    for smi in test_smiles:
        result, reason = is_trichlorobenzene(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")