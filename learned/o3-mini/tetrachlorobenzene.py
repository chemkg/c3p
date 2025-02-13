"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition: Any member of the class of chlorobenzenes carrying four chloro groups 
at unspecified positions. This classifier looks for ANY six‐membered aromatic ring 
in the molecule (even if fused or linked) that has exactly four chlorine substituents 
attached (neighbors not part of the ring). Other substituents are allowed.
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if the given molecule (SMILES string) qualifies as a tetrachlorobenzene.
    The approach is to find any aromatic benzene ring (6-membered, all carbons, aromatic)
    in the molecule and count the number of chlorine (Cl) substituents attached to it 
    (only counting groups on atoms not in the ring). If exactly 4 chlorines are found,
    the molecule is classified as tetrachlorobenzene.

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
        # Ensure the molecule is sanitized and kekulized (aromaticity perception, etc.)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization error: " + str(e)

    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    # List to hold candidate benzene rings (set of atom indices)
    candidate_rings = []
    for ring in ring_info:
        if len(ring) != 6:
            continue
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check that the atom is a carbon and is aromatic
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if is_benzene:
            candidate_rings.append(set(ring))

    if not candidate_rings:
        return False, "No benzene ring (6-membered aromatic carbon cycle) found"

    # For each candidate benzene ring, count the number of Cl substituents.
    for candidate in candidate_rings:
        chlorine_count = 0
        # For every atom in the candidate ring...
        for idx in candidate:
            atom = mol.GetAtomWithIdx(idx)
            # Consider only neighbors not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in candidate:
                    continue
                # Count if neighbor atom is Cl (atomic number 17)
                if nbr.GetAtomicNum() == 17:
                    chlorine_count += 1
        # If exactly four chlorines are attached, we consider it a hit.
        if chlorine_count == 4:
            return True, "Found a benzene ring (aromatic 6-membered ring) with exactly 4 chlorine substituents"
    
    return False, "No benzene ring found with exactly 4 chlorine substituents"


# For testing (these test calls may be removed or commented out)
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
        ("Oc1c(Cl)c(Cl)c(O)c(Cl)c1Cl", "tetrachlorohydroquinone"),
        ("Clc1ccc(cc1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2,3,4,4',5-Pentachlorobiphenyl"),
        ("Clc1cc(Cl)c(Cl)cc1Cl", "1,2,4,5-tetrachlorobenzene"),
        ("Clc1ccc(cc1Cl)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2,3,3',4,4',5-Hexachlorobiphenyl"),
        ("Clc1c(Cl)c(Cl)c2C(=O)OC(=O)c2c1Cl", "tetrachlorophthalic anhydride"),
        ("Clc1c(Cl)c(C#N)c(Cl)c(C#N)c1Cl", "chlorothalonil"),
        ("[O-][N+](=O)c1c(Cl)c(Cl)cc(Cl)c1Cl", "tecnazene"),
    ]
    for smi, name in test_examples:
        result, reason = is_tetrachlorobenzene(smi)
        print(f"{name}: {result} ({reason})")