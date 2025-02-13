"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Definition: Any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
A candidate is identified by finding an isolated (non‐fused) benzene ring (six aromatic carbons) that has exactly three 
terminal chlorine substituents attached directly. In order to reduce false positives (molecules that contain an isolated benzene ring with three Cl’s but where one or more non‐Cl substituents are large and likely indicate that the benzene is only part of a larger scaffold), 
we further inspect non‐chlorine, non‐aromatic substituents: if such a substituent extends beyond a “small” fragment (here, >2 heavy atoms), we disqualify the candidate.
"""

from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule qualifies as a trichlorobenzene based on its SMILES string.
    
    The algorithm:
      1. Parse and sanitize the molecule.
      2. Find six‐membered rings whose atoms are all aromatic carbons (a candidate benzene).
      3. Exclude rings that are fused (i.e. share 2+ atoms with any other ring).
      4. For each candidate ring, examine each ring atom’s neighbors that are not part of the ring.
         • Count a neighbor as a valid chloro substituent if it is a terminal Cl (atomic number 17, degree 1).
         • For any other substituent that is not part of another aromatic ring, compute the “substituent size” 
           (number of heavy atoms reached in the fragment attached through that bond, obtained by a simple DFS). 
           If a non-chloro substituent has more than 2 heavy atoms, we mark this ring as being substituted “complexly.”
      5. If a candidate ring has exactly 3 valid terminal Cl substituents (and no disqualifying complex non‐Cl substituents from non‐aromatic fragments), then we return True.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if an isolated benzene ring with exactly three terminal chlorine substituents is found 
              (and other substituents attached via non‐aromatic bonds are “small”), False otherwise.
        str: Explanation of the classification decision.
    """
    # Helper: recursively count heavy atoms (atomic number > 1) in a substituent fragment starting from a given atom,
    # while excluding coming back to the ring atom.
    def substituent_size(atom, exclude_idx, visited):
        count = 0
        # Only count heavy atoms (ignore hydrogens)
        if atom.GetAtomicNum() > 1:
            count += 1
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == exclude_idx or nbr.GetIdx() in visited:
                continue
            count += substituent_size(nbr, exclude_idx, visited)
        return count

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Loop over rings looking for six-membered aromatic rings (candidate benzene rings)
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        
        # Check every atom in the ring: must be carbon and aromatic 
        is_candidate = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_candidate = False
                break
        if not is_candidate:
            continue
        
        # Exclude fused rings: if this ring shares 2 or more atoms with any other ring, skip it.
        fused = False
        for other_ring in atom_rings:
            if other_ring == ring:
                continue
            if len(set(ring).intersection(other_ring)) >= 2:
                fused = True
                break
        if fused:
            continue
        
        # Now for this candidate ring, count the valid terminal Cl substituents.
        cl_count = 0
        complex_substituent_flag = False
        # Loop over each atom in candidate ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at each neighbor that is not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # If the neighbor is part of an aromatic ring, then assume this bond is the link to another aromatic unit 
                # (and ignore it from substituent complexity checking).
                if nbr.GetIsAromatic():
                    continue
                # If neighbor is chlorine and terminal (degree=1), count as a valid chloro substituent.
                if nbr.GetAtomicNum() == 17 and nbr.GetDegree() == 1:
                    cl_count += 1
                else:
                    # For any other substituent, compute its fragment size.
                    frag_size = substituent_size(nbr, atom.GetIdx(), set())
                    # Allow only very small fragments (threshold = 2 heavy atoms)
                    if frag_size > 2:
                        complex_substituent_flag = True
        # If candidate ring has exactly 3 terminal Cl substituents and no disqualifying substituents, classify as trichlorobenzene.
        if cl_count == 3 and not complex_substituent_flag:
            return True, "Found an isolated benzene ring with exactly three terminal chloro substituents and only small additional substituents"
    
    return False, "No isolated benzene ring with exactly three terminal chloro substituents (with acceptable small substituents) found"


# (Optional) For testing purposes you can use the following:
if __name__ == '__main__':
    test_smiles = [
        "Clc1cccc(Cl)c1Cl",  # 1,2,3-trichlorobenzene (should be True)
        "ClC1=C(OCCN(C(=O)N)CCC)C(Cl)=CC(=C1)Cl"  # BTS44595; expected to be False
    ]
    for smi in test_smiles:
        result, reason = is_trichlorobenzene(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")