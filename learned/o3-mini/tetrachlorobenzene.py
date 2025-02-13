"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition: Any member of the class of chlorobenzenes carrying four chloro groups at unspecified positions.
This improved classifier first identifies candidate benzene rings (6-membered aromatic carbocycles),
then eliminates those that are fused with another benzene ring (sharing ≥2 atoms) and finally checks
each candidate’s substituents. For each ring atom the only allowed substituents are:
   – chlorine (Cl) (which we count),
   – an oxygen that looks like a hydroxyl (i.e. having at least one implicit hydrogen), or 
   – a carbon that is the linkage to another isolated benzene ring (e.g. in biphenyl systems).
If a candidate benzene ring carries exactly 4 chlorine substituents, it is accepted.
Note: This is an heuristic approach and may not be perfect.
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if the given molecule (SMILES string) qualifies as a tetrachlorobenzene.
    The idea is to find at least one isolated benzene ring (6-membered aromatic carbocycle)
    for which the only substituents (non-ring neighbors) are either:
      - Cl atoms (counted toward the four required),
      - O atoms with at least one implicit hydrogen (i.e. hydroxyl groups),
      - or a carbon that is part of another isolated benzene ring (as in biphenyls).
    The candidate ring must have exactly 4 chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a tetrachlorobenzene, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization error: " + str(e)
    
    # Get all rings (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    # Build candidate benzene rings: 6-membered rings where every atom is aromatic carbon.
    benzene_candidates = []
    for ring in ring_info:
        if len(ring) != 6:
            continue
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # check atomic number 6 and aromatic
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if is_benzene:
            benzene_candidates.append(set(ring))
    
    if not benzene_candidates:
        return False, "No 6-membered aromatic ring found"

    # Exclude candidates that are fused with another benzene candidate.
    # (We consider a ring 'fused' if it shares 2 or more atoms with a different benzene candidate.)
    isolated_candidates = []
    for i, ring1 in enumerate(benzene_candidates):
        fused = False
        for j, ring2 in enumerate(benzene_candidates):
            if i == j:
                continue
            if len(ring1.intersection(ring2)) >= 2:
                fused = True
                break
        if not fused:
            isolated_candidates.append(ring1)
    if not isolated_candidates:
        return False, "No isolated benzene ring found (non-fused benzene required)"
    
    # Utility: check if a neighbor atom is in any isolated benzene ring (other than candidate)
    def is_linked_to_benzene(neighbor, candidate_ring):
        # Return True if this neighbor (a carbon) is a member of any isolated candidate ring
        # different from candidate_ring.
        if neighbor.GetAtomicNum() != 6:
            return False
        for ring in isolated_candidates:
            # if neighbor is in a candidate ring and that ring is not the candidate ring being considered
            if neighbor.GetIdx() in ring and ring != candidate_ring:
                return True
        return False

    # Now loop over each isolated (non-fused) benzene candidate and check its substituents.
    for candidate in isolated_candidates:
        chlorine_count = 0
        candidate_valid = True
        # For every atom in the candidate ring...
        for idx in candidate:
            atom = mol.GetAtomWithIdx(idx)
            # Get all neighbors of this ring atom that are not in the ring candidate.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in candidate:
                    continue
                # Allow chloride: count it.
                if nbr.GetAtomicNum() == 17:
                    chlorine_count += 1
                # Allow oxygen if it appears to be -OH (has at least one implicit H)
                elif nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    continue
                # Allow a carbon if it is part of another isolated benzene ring (e.g. biphenyl connection)
                elif nbr.GetAtomicNum() == 6 and is_linked_to_benzene(nbr, candidate):
                    continue
                else:
                    # found an extra substituent (e.g. an acyl group) that disqualifies this candidate
                    candidate_valid = False
                    break
            if not candidate_valid:
                break
        # Accept candidate if it is valid and has exactly 4 chlorine substituents.
        if candidate_valid and chlorine_count == 4:
            return True, "Found an isolated benzene ring with exactly 4 chlorine substituents"
    
    return False, "No isolated benzene ring found with exactly 4 chlorine substituents"


# For testing (you can remove or comment out these test calls)
if __name__ == "__main__":
    # Some examples:
    test_examples = [
        ("Oc1cccc(c1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2',3',4',5'-tetrachlorobiphenyl-3-ol"),
        ("C1(O)=C(C(=C(C=C1Cl)Cl)Cl)Cl", "2,3,4,6-tetrachlorophenol"),
        ("Clc1cccc(c1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2,3,3',4,5-pentachlorobiphenyl"),
        ("Clc1cc(Cl)c(Cl)c(Cl)c1", "1,2,3,5-tetrachlorobenzene"),
        ("COC(=O)c1c(Cl)c(Cl)c(C(=O)OC)c(Cl)c1Cl", "Dacthal (should not be tetrachlorobenzene)"),
        ("Clc1cc2c(oc3c(Cl)c(Cl)c(Cl)c(Cl)c23)c(Cl)c1Cl", "Heptachlorodibenzofuran (false positive in prior method)")
    ]
    for smi, name in test_examples:
        result, reason = is_tetrachlorobenzene(smi)
        print(f"{name}: {result} ({reason})")