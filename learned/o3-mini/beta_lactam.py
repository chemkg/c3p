"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: Beta-lactam
A beta-lactam is defined as a lactam in which the amide bond is contained 
within a four-membered ring. In this context the ring should contain exactly
one nitrogen and three carbons, and one of those carbons (the carbonyl carbon)
must be bonded to a nitrogen (inside the ring) while also bearing an exocyclic
carbonyl (C=O) group. This improved approach is aimed at reducing false positives.
"""

from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    A beta-lactam must contain a four-membered ring composed exactly of one 
    nitrogen and three carbons. In this ring, at least one of the carbon atoms 
    (the carbonyl carbon) must:
       - be directly bonded to a nitrogen atom (within the ring) and 
       - have an external double bond to an oxygen (i.e. an exocyclic carbonyl group)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as beta-lactam, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES and check if valid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If the molecule consists of multiple fragments, check each one.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if not frags:
        frags = [mol]
    
    # Iterate over each fragment.
    for frag in frags:
        # It is often useful to try kekulization to make double bonds explicit.
        try:
            Chem.Kekulize(frag, clearAromaticFlags=True)
        except Exception:
            pass
        
        ring_info = frag.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        
        # Loop over each detected ring.
        for ring in atom_rings:
            if len(ring) != 4:
                continue  # only interested in 4-membered rings
            # Get atoms in ring.
            atoms_in_ring = [frag.GetAtomWithIdx(idx) for idx in ring]
            # Count nitrogen and carbon atoms.
            n_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            c_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 6)
            if n_count != 1 or c_count != 3:
                continue  # does not match beta-lactam ring composition
            
            # Now look for a carbon atom in the ring that meets our criteria:
            # (a) It is bound inside the ring to a nitrogen atom.
            # (b) It has a double bond to an oxygen that is not part of the ring.
            for idx in ring:
                atom = frag.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:  # consider only carbon atoms
                    continue
                
                # Check that at least one neighbor within the ring is the nitrogen.
                neighbors = atom.GetNeighbors()
                has_ring_nitrogen = False
                for nbr in neighbors:
                    if nbr.GetIdx() in ring and nbr.GetAtomicNum() == 7:
                        has_ring_nitrogen = True
                        break
                if not has_ring_nitrogen:
                    continue  # this carbon is not connected to a ring nitrogen
                
                # Now look for a double bond to oxygen on this carbon
                carbonyl_found = False
                for bond in atom.GetBonds():
                    # We require a double bond.
                    if bond.GetBondType() != Chem.BondType.DOUBLE:
                        continue
                    # Get the neighbor bonded by the double bond.
                    nbr = bond.GetOtherAtom(atom)
                    # Check that it is oxygen.
                    if nbr.GetAtomicNum() != 8:
                        continue
                    # and ensure that the oxygen atom is outside the ring.
                    if nbr.GetIdx() in ring:
                        continue
                    # If all conditions hold, then we consider this carbon as the beta-lactam carbonyl.
                    carbonyl_found = True
                    break
                
                # If we found such a carbon atom, then we have a valid beta-lactam ring.
                if carbonyl_found:
                    return True, ("Found beta-lactam ring: a 4-membered ring with exactly one nitrogen and three carbons, "
                                  "in which a carbon (bonded to an in‐ring nitrogen) bears an exocyclic carbonyl group")
                    
    # If no valid beta-lactam ring was discovered:
    return False, "No beta-lactam ring (4-membered lactam with one N and three C having an external carbonyl) found"

# Example usage:
if __name__ == "__main__":
    # Test with oxacillin(1-) SMILES (should be True)
    test_smiles = "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)c1c(C)onc1-c1ccccc1)C([O-])=O"
    result, reason = is_beta_lactam(test_smiles)
    print("Result:", result)
    print("Reason:", reason)