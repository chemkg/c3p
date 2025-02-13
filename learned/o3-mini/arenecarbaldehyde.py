"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
#!/usr/bin/env python
"""
Classifies: Arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is directly attached to an arene.
An arenecarbaldehyde must have a [CX3H1](=O) group that is bonded to an aromatic carbon
which is part of an aromatic ring containing only carbon atoms.
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as an aldehyde in which the carbonyl (C=O) group is directly
    attached to an aromatic carbon that is part of an arene (i.e. a ring composed exclusively of carbon atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an arenecarbaldehyde, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS to match an aldehyde group [CX3H1](=O) bonded to an aromatic atom ([c]).
    # This pattern will return a match where the first atom is the aldehyde carbon and the second is the aromatic neighbor.
    aldehyde_arene_smarts = "[CX3H1](=O)[c]"
    pattern = Chem.MolFromSmarts(aldehyde_arene_smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern"

    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No aldehyde group directly attached to an aromatic atom was found"

    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()

    # For each aldehyde match, further check that the aromatic neighbor is indeed a carbon atom and
    # that it belongs to at least one aromatic ring that is composed solely of carbon atoms.
    for match in matches:
        # match[0] is the aldehyde carbon; match[1] is its aromatic neighbor.
        aromatic_neighbor = mol.GetAtomWithIdx(match[1])
        # Ensure the neighbor is a carbon; if not, skip this match.
        if aromatic_neighbor.GetAtomicNum() != 6:
            continue

        # Find all rings that contain the aromatic_neighbor.
        neighbor_idx = aromatic_neighbor.GetIdx()
        rings = ring_info.AtomRings()
        pure_arene_found = False
        for ring in rings:
            if neighbor_idx in ring:
                # Check that every aromatic atom in the ring is a carbon.
                all_carbons = True
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetIsAromatic() and atom.GetAtomicNum() != 6:
                        all_carbons = False
                        break
                if all_carbons:
                    pure_arene_found = True
                    break
        if pure_arene_found:
            return True, "Found an aldehyde group with the carbonyl carbon attached to a pure arene ring (all aromatic atoms are carbon)."

    return False, "Aldehyde group found, but the attached aromatic atom is not part of a pure arene ring."

# Uncomment below for sample testing:
# test_smiles = [
#     "[H]C(=O)c1cccc2ccccc12",  # 1-naphthaldehyde (should be True)
#     "Cc1cccc(C=O)c1",          # m-tolualdehyde (should be True)
#     "CN1C2=CC=CC=C2C(=C1SC3=CC=C(C=C3)Cl)C=O",  # false positive in previous version (should be False)
#     "[O-][N+](=O)c1ccc(C=O)cc1",  # 4-nitrobenzaldehyde (should be True)
# ]
# for smi in test_smiles:
#     result, reason = is_arenecarbaldehyde(smi)
#     print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")