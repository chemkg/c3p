"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition: Any member of the class of chlorobenzenes carrying four chloro groups
at unspecified positions. This improved classifier looks for a benzene ring (6-membered aromatic carbon cycle)
that has exactly four chlorine (Cl) substituents attached (neighbors not in the ring).
To avoid false positives from highly fused polycyclic aromatic systems, it favors rings that are not fused 
(with another benzene ring, i.e. sharing 2 or more atoms) unless none are found.
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines whether the given molecule (SMILES string) qualifies as a tetrachlorobenzene.
    The classifier first looks for benzene rings (aromatic 6-membered rings composed only of carbons).
    It then counts the number of chlorine substituents (neighbors not belonging to the ring).
    In order to reduce false positives from polycyclic (fused) systems, we attempt to filter out 
    any candidate ring that is fused with another benzene ring (i.e. shares 2 or more atoms with it),
    as long as at least one isolated benzene is present.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a benzene ring with exactly 4 chlorine substituents is found, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization error: " + str(e)
    
    # Retrieve ring information
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []
    for ring in ring_info:
        # Only consider rings of exactly 6 atoms.
        if len(ring) != 6:
            continue
        # Check if every atom in the ring is carbon and aromatic (benzene-like)
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if is_benzene:
            candidate_rings.append(set(ring))
    
    if not candidate_rings:
        return False, "No 6‐membered aromatic carbon ring (benzene) found."
    
    # Try to filter out rings that are fused with another benzene ring.
    # Two benzene rings that share 2 or more atoms are “fused”
    isolated_rings = []
    for i, ring in enumerate(candidate_rings):
        fused = False
        for j, other in enumerate(candidate_rings):
            if i == j:
                continue
            if len(ring.intersection(other)) >= 2:
                fused = True
                break
        if not fused:
            isolated_rings.append(ring)
    
    # If at least one isolated benzene ring is available, use that list; otherwise, fall back.
    rings_to_check = isolated_rings if isolated_rings else candidate_rings

    # For each candidate ring, count Cl substituents attached to its atoms (neighbors not in ring).
    for ring in rings_to_check:
        chlorine_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms that are part of the ring
                if nbr.GetAtomicNum() == 17:  # chlorine atomic number is 17
                    chlorine_count += 1
        if chlorine_count == 4:
            return True, "Found a benzene ring with exactly 4 chlorine substituents."
    
    return False, "No benzene ring found with exactly 4 chlorine substituents on its periphery."

# For manual testing (remove or comment out in production)
if __name__ == "__main__":
    test_examples = [
        ("Oc1cccc(c1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2',3',4',5'-tetrachlorobiphenyl-3-ol"),
        ("C1(O)=C(C(=C(C=C1Cl)Cl)Cl)Cl", "2,3,4,6-tetrachlorophenol"),
        ("Clc1cccc(c1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2,3,3',4,5-pentachlorobiphenyl"),
        ("Clc1cc(Cl)c(Cl)c(Cl)c1", "1,2,3,5-tetrachlorobenzene"),
        ("C1(O)=C(C(=C(C(=C1)Cl)Cl)Cl)Cl", "2,3,4,5-tetrachlorophenol"),
        ("Oc1c(Cl)c(Cl)cc(Cl)c1Cl", "2,3,5,6-tetrachlorophenol"),
        ("Clc1ccc(Cl)c(Cl)c1Cl", "1,2,3,4-tetrachlorobenzene"),
        ("Clc1cc(-c2ccccc2)c(Cl)c(Cl)c1Cl", "2,3,4,5-tetrachlorobiphenyl"),
        ("Oc1c(O)c(Cl)c(Cl)c(Cl)c1Cl", "tetrachlorocatechol"),
        ("Clc1c(Cl)c(C#N)c(Cl)c(C#N)c1Cl", "chlorothalonil"),
        ("[O-][N+](=O)c1c(Cl)c(Cl)cc(Cl)c1Cl", "tecnazene"),
    ]
    for smi, name in test_examples:
        result, reason = is_tetrachlorobenzene(smi)
        print(f"{name}: {result} ({reason})")