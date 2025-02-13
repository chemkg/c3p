"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: Beta-lactam
A beta-lactam is defined as a lactam (cyclic amide) in which the amide bond is contained 
within a four-membered ring. In practice, beta-lactams have a 4-membered ring made up of exactly 
one nitrogen and three carbon atoms, where one of these carbons is part of a carbonyl group 
(double-bonded to an oxygen that lies outside the ring).
"""

from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam must contain a four-membered ring composed exactly of one nitrogen and 
    three carbons, with one of the carbon atoms bearing a carbonyl group (double bond O that is external to the ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as beta-lactam, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If there are multiple fragments (e.g. mixtures), iterate over each one.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if not frags:
        frags = [mol]
    
    for frag in frags:
        try:
            # Kekulize to ensure double bonds (especially C=O) are explicit.
            Chem.Kekulize(frag, clearAromaticFlags=True)
        except Exception:
            # Continue even if kekulization fails.
            pass

        ring_info = frag.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        # Iterate over each ring present in the fragment.
        for ring in atom_rings:
            # Only consider rings of exactly 4 atoms.
            if len(ring) != 4:
                continue

            # Extract atoms belonging to the ring.
            atoms_in_ring = [frag.GetAtomWithIdx(idx) for idx in ring]
            
            # Count number of nitrogen and carbon atoms in the ring.
            n_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            c_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 6)
            
            # Beta-lactam rings should contain exactly one nitrogen and three carbons.
            if n_count != 1 or c_count != 3:
                continue

            # Check if one of the carbon atoms is carbonyl-like.
            carbonyl_found = False
            for atom in atoms_in_ring:
                if atom.GetAtomicNum() != 6:
                    continue
                # Look at each bond for a double bond to oxygen.
                for bond in atom.GetBonds():
                    # Verify bond is a double bond.
                    if bond.GetBondType() != Chem.BondType.DOUBLE:
                        continue
                    # Get the neighboring atom.
                    neighbor = bond.GetOtherAtom(atom)
                    # Check for oxygen.
                    if neighbor.GetAtomicNum() != 8:
                        continue
                    # Make sure that the oxygen is not in the current ring.
                    if neighbor.GetIdx() in ring:
                        continue
                    carbonyl_found = True
                    break
                if carbonyl_found:
                    break
            
            if carbonyl_found:
                return True, ("Found beta-lactam ring: a 4-membered ring with exactly one nitrogen and three carbons, "
                              "with one carbon attached via a double bond to an oxygen external to the ring")
    
    # No qualifying beta-lactam ring found across all fragments.
    return False, "No beta-lactam ring (4-membered lactam with one N and three C with external carbonyl) found"

# Example usage:
if __name__ == "__main__":
    # Test with oxacillin(1-) SMILES (should be True)
    test_smiles = "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)c1c(C)onc1-c1ccccc1)C([O-])=O"
    result, reason = is_beta_lactam(test_smiles)
    print("Result:", result, "\nReason:", reason)