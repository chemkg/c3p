"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: Azole
Definition: Any heteroarene with a five-membered aromatic ring that contains at least one nitrogen.
Note: Although the original definition said "monocyclic," many azole drug molecules
contain a fused system. Here we look for any five-membered aromatic ring containing nitrogen.
We also try to avoid classifying peptides (which often contain histidine azole rings)
by a simple heuristic.
"""

from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as a molecule containing a five‐membered aromatic ring with at least one nitrogen.
    In our approach we allow for fused ring systems.
    Also, if the molecule appears to be a peptide (by having two or more amide bonds and a low
    heavy-atom count) we reject classifying it as an azole even if it contains an azole ring.

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        bool: True if the molecule contains a qualifying azole ring, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize the molecule to update aromaticity and ring info.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {str(e)}"

    # Get all ring atom index sets in the molecule.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in the molecule"

    # Look for any 5-membered ring that is aromatic and contains at least one nitrogen.
    azole_found = False
    for ring in rings:
        if len(ring) != 5:
            continue

        # Check that every atom in the ring is marked aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # Check if at least one nitrogen is present in this ring.
        if not any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
            continue

        # We have found a five-membered aromatic ring with nitrogen.
        azole_found = True
        break

    if not azole_found:
        return False, "No qualifying five-membered aromatic ring containing nitrogen found"

    # Heuristic to reject peptides:
    # Many peptides contain one or more azole (histidine) rings but should not be classified as azoles.
    # We check if the molecule has two or more amide bonds and few heavy atoms.
    peptide_smarts = Chem.MolFromSmarts("[CX3](=O)[NX3]")  # simple amide bond pattern
    amide_matches = mol.GetSubstructMatches(peptide_smarts)
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if len(amide_matches) >= 2 and num_heavy_atoms < 50:
        return False, "Molecule appears to be a peptide with an azole side chain"

    return True, "Found a five-membered aromatic ring containing nitrogen (azole) in the molecule"


# Example usage when running as a script:
if __name__ == "__main__":
    test_smiles = [
        "[N+](C)([C@H](C(=O)[O-])CC=1NC(SC)=NC1)(C)C",  # S-methyl-L-ergothioneine (should be True)
        "O(C(=O)C=1NC=CC1)C",                           # Methyl 1H-pyrrole-2-carboxylate (should be True)
        "C1C=CN=C1",                                   # 3H-pyrrole (should be True even if fused concept not used)
        "O=C(N[C@@H](CC=1NC=NC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)[C@H](O)CC2=CC=CC=C2)",  # Thr-Phe-His (peptide, should be False)
    ]
    
    for smile in test_smiles:
        result, reason = is_azole(smile)
        print(f"SMILES: {smile}\nResult: {result}\nReason: {reason}\n")