"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: Azole
Definition: Any monocyclic heteroarene consisting of a five‐membered ring containing nitrogen.
Azoles can also contain one or more other non‐carbon atoms, such as nitrogen, sulfur or oxygen.
"""

from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as any monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    The five-membered ring must be aromatic and non-fused (each atom in the ring should belong to only that ring).

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        bool: True if the molecule contains a qualifying azole ring, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize the molecule to ensure aromaticity and ring info are up to date.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {str(e)}"

    # Get all ring atom indices as tuples
    ring_atom_tuples = mol.GetRingInfo().AtomRings()
    if not ring_atom_tuples:
        return False, "No rings found in the molecule"

    # For each ring, check if it qualifies as a five-membered azole.
    for ring in ring_atom_tuples:
        if len(ring) != 5:
            continue

        # Check that every atom in this ring is aromatic.
        all_aromatic = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic():
                all_aromatic = False
                break
        if not all_aromatic:
            continue

        # Check that the ring is non-fused.
        # In other words, for each atom in the ring, it should appear only in this ring.
        ring_isolated = True
        for idx in ring:
            count = sum(1 for r in ring_atom_tuples if idx in r)
            if count > 1:
                ring_isolated = False
                break
        if not ring_isolated:
            continue

        # Check for at least one nitrogen atom in the ring.
        has_nitrogen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring)
        if not has_nitrogen:
            continue

        # Found a qualifying five-membered aromatic non-fused ring with nitrogen.
        return True, "Found five-membered aromatic non-fused ring containing nitrogen (azole) in the molecule"

    return False, "No five-membered monocyclic aromatic ring containing nitrogen (azole) found"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "[N+](C)([C@H](C(=O)[O-])CC=1NC(SC)=NC1)(C)C",  # S-methyl-L-ergothioneine (true positive)
        "C1C=CN=C1",  # 3H-pyrrole (should be true)
        "O=C1C2=C(O)C3=C(O)C=4C5=C6N(OC(C6=CC=C5)=O)CC4C=C3O[C@@]2(C(=O)OC)[C@@H](O)CC1"  # complex false positive candidate
    ]
    
    for smile in test_smiles:
        result, reason = is_azole(smile)
        print(f"SMILES: {smile}\nResult: {result}\nReason: {reason}\n")