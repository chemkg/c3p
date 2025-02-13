"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: Compounds containing a 1,2,4-triazine skeleton
Definition: Any compound with a 1,2,4-triazine skeleton, in which nitrogen atoms replace carbon 
at positions 1, 2 and 4 of the core benzene ring structure.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine skeleton.
    A valid 1,2,4-triazine is defined here as a six-membered aromatic ring that contains
    exactly three nitrogen atoms and three carbon atoms. In addition, there must exist a rotation 
    (or its reverse) of the ring order where the pattern is: [N, N, C, N, C, C],
    corresponding to nitrogen atoms at positions 1, 2, and 4 (with positions 1 and 2 adjacent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a 1,2,4-triazine skeleton, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings detected in molecule"

    # Define the target pattern for a 1,2,4-triazine ring.
    # When laid out linearly (with cyclic rotation) we want: N, N, C, N, C, C.
    target_forward = ["N", "N", "C", "N", "C", "C"]
    target_reverse = list(reversed(target_forward))
    
    # Iterate over all rings
    for ring in ring_info:
        if len(ring) != 6:
            continue  # we need a 6-membered ring

        # Get the atoms in the ring in the given order (as provided by RDKit)
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check that all atoms in the ring are aromatic.
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue

        # Check that the ring has only C and N atoms.
        symbols = [atom.GetSymbol() for atom in ring_atoms]
        if any(sym not in ["C", "N"] for sym in symbols):
            continue

        # Count the number of nitrogen atoms in the ring.
        if symbols.count("N") != 3 or symbols.count("C") != 3:
            continue

        # Try all rotations (and the reverse order) to see if the pattern matches.
        n = len(symbols)
        matched = False
        for i in range(n):
            # rotation of the symbols list:
            rotation = symbols[i:] + symbols[:i]
            if rotation == target_forward or rotation == target_reverse:
                matched = True
                break

        if matched:
            return True, "Contains a 1,2,4-triazine skeleton"
    
    # If no ring passed the tests that qualify the 1,2,4-triazine skeleton:
    return False, "No 1,2,4-triazine skeleton found"