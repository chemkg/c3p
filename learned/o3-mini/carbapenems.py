"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: Carbapenems – beta-lactam antibiotics with a carbapenem skeleton substituted at positions 3, 4, and 6.
This improved heuristic implementation detects a four-membered beta-lactam ring (with one nitrogen and at least one carbonyl group)
fused with a five-membered ring that (a) does not contain sulfur and (b) includes at least one carbon-carbon double bond,
a characteristic of the unsaturated five-membered ring in many carbapenem antibiotics.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    This heuristic implementation works as follows:
      1. Parse the SMILES and get ring information.
      2. Look for candidate 4-membered rings (potential beta-lactam rings) that have:
            - Exactly 4 atoms.
            - Exactly one nitrogen atom.
            - At least one carbon with a double-bonded oxygen (i.e. a carbonyl).
      3. For each candidate 4-membered ring, look for a fused 5-membered ring that:
            - Shares at least 2 atoms with the 4-membered ring.
            - Does not contain any sulfur.
            - Contains at least one carbon-carbon double bond inside the ring.
      4. If such a fused bicyclic system is found, classify the molecule as a carbapenem.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a carbapenem, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings present in the molecule"
    
    carbapenem_core_found = False
    reason = ""
    
    # Iterate over all rings to look for candidate 4-membered rings (beta-lactam rings)
    for ring in ring_info:
        if len(ring) != 4:
            continue
        
        # Check: The 4-ring must contain exactly one nitrogen.
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen != 1:
            continue

        # Check for at least one carbon in the ring having a double-bonded oxygen (carbonyl group)
        has_carbonyl = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # carbon
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:  # oxygen
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            has_carbonyl = True
                            break
            if has_carbonyl:
                break
        if not has_carbonyl:
            continue

        # Now search for a fused 5-membered ring
        for other_ring in ring_info:
            if len(other_ring) != 5:
                continue
            # Check if the candidate 4-membered ring and this 5-membered ring share at least 2 atoms (fused)
            common_atoms = set(ring).intersection(other_ring)
            if len(common_atoms) < 2:
                continue

            # Exclude rings that contain sulfur (atomic num 16). Usually penicillins have sulfur.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 16 for idx in other_ring):
                continue

            # Check that the five-membered ring has at least one carbon-carbon double bond.
            has_c_c_double = False
            # Loop over bonds in the 5-membered ring.
            for i in range(len(other_ring)):
                idx1 = other_ring[i]
                idx2 = other_ring[(i+1) % len(other_ring)]
                bond = mol.GetBondBetweenAtoms(idx1, idx2)
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    atom1 = mol.GetAtomWithIdx(idx1)
                    atom2 = mol.GetAtomWithIdx(idx2)
                    if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                        has_c_c_double = True
                        break
            if not has_c_c_double:
                continue

            # If we reach here, we found a 4-membered ring with proper features fused with a 5-membered ring that has a C=C bond.
            carbapenem_core_found = True
            reason = ("Molecule contains a beta-lactam ring (4-membered; one N and a carbonyl) fused with a "
                      "5-membered unsaturated (C=C present) ring lacking sulfur, which is consistent with a carbapenem skeleton")
            break
        if carbapenem_core_found:
            break

    if not carbapenem_core_found:
        return False, "No fused bicyclic system (beta-lactam ring fused with an unsaturated, sulfur-free 5-membered ring) found"
    
    # Optional: Further checks on substitution at positions 3,4,6 may be added here.
    return True, reason

# Example usage (you may remove or modify these before production):
if __name__ == "__main__":
    # Test a couple of SMILES strings from the provided examples:
    test_smiles = [
        "[H][C@]12CC=CN1C(=O)C2",  # carbapenem (should be True)
        "OC\\C=C1\\O[C@@H]2CC(=O)N2[C@H]1C(O)=O"  # Isoclavulanic acid (should be False)
    ]
    for s in test_smiles:
        result, msg = is_carbapenems(s)
        print(f"SMILES: {s}\nResult: {result}, Reason: {msg}\n")