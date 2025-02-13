"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: Beta-lactam
A beta-lactam is defined as a lactam in which the amide bond is contained 
within a four-membered ring. In our approach the ring must contain exactly
one nitrogen and three carbons. Moreover, one of the carbon atoms in the ring 
(the “beta‐lactam carbonyl”) must:
   - be connected within the ring to a nitrogen,
   - possess an exocyclic double bond to an oxygen atom (i.e. the oxygen is not part 
     of the ring), and
   - have a total degree of 3 (ensuring that it is not decorated with extra substituents).
This extra check is intended to reduce false positive classifications.
"""

from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    A beta-lactam must contain at least one four-membered ring consisting
    of exactly one nitrogen and three carbons. In that ring one carbon (the 
    beta-lactam carbonyl) must:
       (a) be bonded inside the ring to a nitrogen atom and a carbon atom,
       (b) have an exocyclic double bond to an oxygen atom (the oxygen not in 
           the ring), and
       (c) have a total degree (number of bonded neighbors) equal to 3.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as beta-lactam, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Many molecules contain more than one fragment. Process each fragment.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if not frags:
        frags = [mol]
    
    # Loop over each fragment.
    for frag in frags:
        # Attempt kekulization so that bond orders (especially double bonds) are explicit.
        try:
            Chem.Kekulize(frag, clearAromaticFlags=True)
        except Exception:
            pass
            
        ring_info = frag.GetRingInfo()
        # Iterate over all rings.
        for ring in ring_info.AtomRings():
            if len(ring) != 4:
                continue  # we are only interested in 4-membered rings
            
            # Get the atoms in the ring.
            atoms_in_ring = [frag.GetAtomWithIdx(idx) for idx in ring]
            # Count nitrogen and carbon atoms.
            n_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            c_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 6)
            if n_count != 1 or c_count != 3:
                continue  # the ring composition does not match a beta-lactam
            
            # For each carbon in the ring, check if it can be the beta-lactam carbonyl.
            for idx in ring:
                atom = frag.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue  # only consider carbon atoms
                
                # Check that this carbon has a ring neighbor that is nitrogen.
                ring_nitrogen_found = False
                ring_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring]
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring and nbr.GetAtomicNum() == 7:
                        ring_nitrogen_found = True
                        break
                if not ring_nitrogen_found:
                    continue  # this carbon is not directly bonded to the ring nitrogen
                
                # Now check for an exocyclic double bond to an oxygen.
                # Also, verify that the carbon's degree is exactly 3 
                # (i.e. it has no extra substituents beyond the ring and the carbonyl oxygen).
                if atom.GetDegree() != 3:
                    continue
                
                exocyclic_oxygens = []
                for bond in atom.GetBonds():
                    # Skip bonds to atoms that are in the ring.
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetIdx() in ring:
                        continue
                    # We require the bond to be a double bond.
                    if bond.GetBondType() != Chem.BondType.DOUBLE:
                        continue
                    # Check that the neighbor is an oxygen.
                    if nbr.GetAtomicNum() == 8:
                        exocyclic_oxygens.append(nbr)
                # If we found at least one double-bonded exocyclic oxygen, return True.
                if exocyclic_oxygens:
                    return True, ("Found beta-lactam ring: a 4-membered ring with exactly one nitrogen and three carbons, "
                                  "where a carbon (bonded to an in‐ring nitrogen) bears an exocyclic carbonyl group "
                                  "and shows no extra substituents (degree 3)")
    return False, "No beta-lactam ring (4-membered lactam with one N and three C having the correct carbonyl) found"

# Example usage (for testing):
if __name__ == "__main__":
    # Using oxacillin(1-) as a test (should be classified as beta-lactam):
    test_smiles = "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)c1c(C)onc1-c1ccccc1)C([O-])=O"
    result, reason = is_beta_lactam(test_smiles)
    print("Result:", result)
    print("Reason:", reason)