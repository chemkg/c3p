"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: Carbapenems – beta-lactam antibiotics with a carbapenem skeleton substituted at positions 3, 4, and 6.
This heuristic implementation detects a four-membered beta-lactam ring (containing one nitrogen and a carbonyl)
that is fused with a five-membered ring that does not contain sulfur.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    The detection is heuristic: we search for a four-membered beta-lactam ring with one nitrogen and one carbonyl,
    fused (sharing at least 2 atoms) with a five-membered ring that does not contain sulfur (differentiating it from penicillins).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a carbapenem, False otherwise
        str: Reason for the classification decision
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information (as tuples of atom indices for each ring)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings present in the molecule"

    # Initialize a flag for detecting a fused beta-lactam system
    carbapenem_core_found = False

    # First, search for candidate 4-membered rings (potential beta-lactam rings)
    for ring in ring_info:
        if len(ring) != 4:
            continue

        # Count the number of nitrogen atoms in the ring.
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen != 1:
            continue

        # Look for at least one carbon within the ring that bears a double-bonded oxygen (carbonyl)
        has_carbonyl = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # carbon
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            has_carbonyl = True
                            break
            if has_carbonyl:
                break
        if not has_carbonyl:
            continue

        # Now, check if this 4-membered ring is fused with a 5-membered ring.
        # Fusion is assumed if at least 2 atoms are shared with a 5-membered ring.
        fused_with_5 = False
        for other_ring in ring_info:
            if len(other_ring) != 5:
                continue
            # Check if the rings share at least 2 atoms.
            common_atoms = set(ring).intersection(other_ring)
            if len(common_atoms) >= 2:
                # To help distinguish carbapenems from penicillins,
                # ensure that the five-membered ring does not contain sulfur (atomic num 16).
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() != 16 for idx in other_ring):
                    fused_with_5 = True
                    break
        if fused_with_5:
            carbapenem_core_found = True
            break

    if not carbapenem_core_found:
        return False, "No fused beta-lactam (4-membered ring with carbonyl and one N) and five-membered ring (lacking S) found"

    # Optionally, more detailed checks (such as substitution at positions 3,4,6) could be included.
    return True, "Molecule contains a beta-lactam ring fused with a five-membered ring (carbapenem skeleton) likely with appropriate substitution"