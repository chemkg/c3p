"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition: Any member of the class of chlorobenzenes carrying four chloro groups at unspecified positions.
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if the given molecule (SMILES string) is a tetrachlorobenzene.
    
    The molecule qualifies if it contains at least one benzene ring (6-membered aromatic carbocycle)
    that has exactly four chlorine (Cl) substituents (attached to the ring atoms but not in the ring).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a tetrachlorobenzene, False otherwise.
        str: Reason for Classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule has aromaticity perceived.
    Chem.SanitizeMol(mol)

    # Get ring information
    ring_info = mol.GetRingInfo().AtomRings()

    # Iterate over each ring and check if it is an aromatic benzene ring that carries 4 Cl substituents.
    for ring in ring_info:
        # Check if ring size is 6.
        if len(ring) != 6:
            continue

        # Check that all ring atoms are carbon and aromatic.
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue

        # Count chlorine substituents attached to ring atoms (neighbors that are not in the ring).
        chlorine_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring and nbr.GetAtomicNum() == 17:  # 17 is the atomic number of Cl
                    chlorine_count += 1

        if chlorine_count == 4:
            return True, "Found a benzene ring with exactly 4 chlorine substituents"

    return False, "No benzene ring with exactly 4 chlorine substituents found"