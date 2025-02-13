"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies amino sugars:
Definition: Any sugar having one or more alcoholic hydroxy groups replaced by substituted or unsubstituted amino groups.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    Strategy:
      1. Parse the SMILES string.
      2. Check for the presence of a sugar-like ring – for instance a 5- or 6-membered ring that contains at least one oxygen.
      3. Look for an amino substituent (–NH2) attached to a ring carbon.
         This simple heuristic uses a SMARTS substructure pattern [C;R]-[NH2],
         and then skips cases where the amino nitrogen is part of an amide (i.e. directly attached to a carbonyl).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an amino sugar, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use ring info: search for a candidate sugar ring (5- or 6-membered ring containing at least one oxygen)
    ring_info = mol.GetRingInfo()
    candidate_ring_found = False
    atom_indices_in_candidate_ring = set()
    for ring in ring_info.AtomRings():
        if len(ring) in [5, 6]:
            # check if one of the atoms in the ring is oxygen
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            if any(atom.GetAtomicNum() == 8 for atom in ring_atoms):
                candidate_ring_found = True
                atom_indices_in_candidate_ring.update(ring)
    if not candidate_ring_found:
        return False, "No sugar-like ring (5- or 6-membered ring with oxygen) detected"

    # Define a SMARTS pattern for an amino substituent
    # This simple pattern will match an sp3 carbon in a ring attached to an NH2 group.
    amino_pattern = Chem.MolFromSmarts("[C;R]-[NH2]")
    if amino_pattern is None:
        return False, "Error setting up amino substructure pattern"

    matches = mol.GetSubstructMatches(amino_pattern)
    if not matches:
        return False, "No amino substituents (NH2 on a ring carbon) found"
    
    # Check the found matches to ensure:
    # 1. The carbon bearing the amino group belongs to one of the candidate sugar rings.
    # 2. The nitrogen in question is not part of an amide (i.e. it does not neighbor a carbon with a double-bonded oxygen).
    valid_amino = False
    for match in matches:
        # match: tuple (carbonIdx, nitrogenIdx)
        c_idx, n_idx = match
        if c_idx not in atom_indices_in_candidate_ring:
            # Only consider amino groups on a ring carbon that is part of the candidate ring.
            continue
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Check bonds from the nitrogen: if any neighboring carbon is bound to an oxygen by a double bond, then it is likely part of an amide.
        is_amide = False
        for nb in n_atom.GetNeighbors():
            if nb.GetAtomicNum() == 6:
                for bond in nb.GetBonds():
                    # find a double bond to oxygen
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other = bond.GetOtherAtom(nb)
                        if other.GetAtomicNum() == 8:
                            is_amide = True
                            break
                if is_amide:
                    break
        if is_amide:
            continue
        # If we reach here, we have found a valid amino substituent on a sugar ring.
        valid_amino = True
        break

    if not valid_amino:
        return False, "Amino substituents found but none on a sugar ring or valid (non-amide) amino group"
    
    # Optionally, one can add further checks such as molecular weight or number of hydroxyl groups
    # to avoid false positives.
    
    return True, "Molecule contains a sugar-like ring with at least one amino substituent replacing an -OH group"

# Below is one example usage (uncomment for testing):
# example_smiles = "N[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O"  # 3-amino-3-deoxy-D-glucopyranose
# print(is_amino_sugar(example_smiles))