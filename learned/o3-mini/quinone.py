"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: Quinones, defined as compounds having a fully conjugated cyclic dione structure.
For example, benzoquinones and their polycyclic or heterocyclic analogues. The heuristic
checks for a ring (of size ≥5) that contains exactly two carbonyl groups embedded in a conjugated system.
"""

from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    The function searches each ring in the molecule for exactly two carbonyl groups
    (C(=O)) that are embedded in a conjugated system. In addition, each carbonyl in the ring
    must be connected (within the ring) to two sp2-hybridized (or aromatic) atoms.
    Also, we require that at least 2/3 of the atoms in the ring are sp2 (or aromatic) so that
    the ring can be considered “fully conjugated”.
    
    This heuristic may yield false positives or negatives but is tuned based on a set of outcomes.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as quinone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the molecule and check validity.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Do NOT clear aromatic flags so that a few rings still carry their original aromatic info.
    # However, we will try to assign hybridization.
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # Continue if sanitization fails (although classification may then be uncertain).
        pass

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Helper: check if an atom is a carbonyl carbon 
    # (i.e. a carbon connected to at least one oxygen via a double bond)
    def is_carbonyl(atom):
        if atom.GetAtomicNum() != 6:
            return False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    return True
        return False

    # For each ring, check our quinone criteria.
    for ring in atom_rings:
        if len(ring) < 5:
            # small rings unlikely to be the conjugated quinone core
            continue

        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Count carbonyl atoms in the ring.
        carbonyl_atoms = [atom for atom in ring_atoms if is_carbonyl(atom)]
        if len(carbonyl_atoms) != 2:
            continue

        # For each carbonyl atom, check that both of its neighbors in the ring are sp2-hybridized or aromatic.
        # (This helps ensure that the carbonyl is embedded in a conjugated system.)
        meets_embedded = True
        for atom in carbonyl_atoms:
            # Find neighbors that are in the same ring.
            n_in_ring = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring]
            # Count how many of these have sp2 hybridization or are aromatic.
            count_conjugated = sum(1 for nbr in n_in_ring if (nbr.GetIsAromatic() or nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP2))
            if count_conjugated < 2:
                meets_embedded = False
                break
        if not meets_embedded:
            continue

        # Check overall ring conjugation: require at least 2/3 of atoms are sp2-hybridized (or aromatic)
        conjugated_count = sum(1 for atom in ring_atoms if (atom.GetIsAromatic() or atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2))
        if conjugated_count < (2/3) * len(ring_atoms):
            continue

        # If the ring has passed the above criteria, we classify as quinone.
        return True, "Found a ring with two carbonyl groups embedded in a fully conjugated cyclic system"
         
    # If no ring qualifies, return false with explanation.
    return False, "No fully conjugated cyclic dione (quinone) structure was found"


# For testing purposes (optional)
if __name__ == "__main__":
    # List a few test SMILES (both quinones and non-quinones)
    test_smiles = [
        "COC1=C(OC)C(=O)C(C)=CC1=O",  # ubiquinone-0 (expected True)
        "O=C1C=C(OC)[C@](O)(C)C([C@@]1(OC(=O)C(CC)C)C)=O",  # Phomaligol A (expected True per example – our new scheme might capture it)
        "C1(=CC=C2NC(C(C2=C1)=O)=O)Br",  # 5-Bromoisatin (false positive previously – we hope to reject)
    ]
    for s in test_smiles:
        result, reason = is_quinone(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*40}")