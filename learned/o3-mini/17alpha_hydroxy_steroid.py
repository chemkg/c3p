"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 17alpha-hydroxy steroid
Definition: The alpha-stereoisomer of 17-hydroxy steroid.
Heuristic (improved):
  1. Require a fused steroid nucleus: exactly one 5-membered ring (D-ring)
     and at least 3 six-membered rings (A, B, C). We then build a steroid_core
     as the union of the D-ring and those six-membered rings that share atoms with it.
  2. In the D-ring, look for a candidate carbon (atomic #6) that is:
       - Substituted by an –OH (an oxygen neighbor bearing at least one hydrogen)
       - Chiral (has a defined chiral tag)
       - Has three heavy-atom (non-H) neighbors (besides the hydroxyl oxygen, there
         should be three carbon neighbors as expected for C17)
       - Has at least one heavy neighbor that is not in the steroid_core (this should
         be the exocyclic side chain at C17)
  3. If exactly one candidate is found, report success (including any CIP assignment),
     otherwise report failure.
Note: This heuristic remains simplified, since the full stereochemical assignment
      of steroids can be complex.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule belongs to the 17alpha-hydroxy steroid class based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 17alpha-hydroxy steroid, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Force assignment of stereochemistry (including CIP codes) if possible.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Get ring information. We require the presence of rings.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings detected – not a steroid"
    
    # Identify the unique 5-membered ring (our D-ring candidate).
    five_membered = [ring for ring in rings if len(ring) == 5]
    if len(five_membered) != 1:
        return False, f"Expected exactly 1 five-membered ring (steroid D-ring), but found {len(five_membered)}"
    
    # Check that there are at least 3 six-membered rings (steroid A, B, C).
    six_membered = [ring for ring in rings if len(ring) == 6]
    if len(six_membered) < 3:
        return False, f"Expected at least 3 six-membered rings (steroid A, B, C), but found {len(six_membered)}"
    
    # Build the fused steroid core: start with the set of atoms in the D-ring,
    # and add any six-membered ring that shares an atom with the D-ring.
    d_ring = five_membered[0]
    steroid_core = set(d_ring)
    for ring in six_membered:
        if set(ring).intersection(steroid_core):
            steroid_core.update(ring)
    # (For our purposes, steroid_core is expected to contain most atoms of the fused steroid nucleus.)
    
    candidate_found = None
    candidate_reason = ""
    
    # Iterate over atoms in the D-ring looking for a candidate C17.
    for atom_idx in d_ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Expect a carbon atom in the D-ring.
        if atom.GetAtomicNum() != 6:
            continue
            
        # Look for an -OH group: an oxygen neighbor that itself has at least one hydrogen.
        oh_found = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                oh_found = True
                break
        if not oh_found:
            continue
        
        # Check that the atom is chiral.
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            candidate_reason = "Hydroxyl-bearing carbon in D-ring is not stereochemically defined"
            continue
        
        # Count heavy-atom neighbors (ignore hydrogens and the oxygen from the OH group separately).
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # For a typical C17, expect three carbon neighbors (one from the fused core and one exocyclic).
        if len(heavy_neighbors) != 4:  
            # Typically, C17 in a steroid is bonded to: two ring carbons, one oxygen (the OH), and one exocyclic carbon.
            # Thus heavy neighbor count including oxygen is 4.
            continue
        
        # Now check that at least one heavy neighbor (apart from the oxygen) lies outside the steroid_core.
        exocyclic = False
        for nbr in heavy_neighbors:
            # Skip oxygen neighbor (the -OH)
            if nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetIdx() not in steroid_core:
                exocyclic = True
                break
        if not exocyclic:
            candidate_reason = "Hydroxyl-bearing carbon in D-ring does not show an exocyclic substituent"
            continue
        
        # If we have reached here, candidate meets our criteria.
        # (Optionally report the CIP code if it exists.)
        cip = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else "unknown"
        candidate_found = (atom_idx, cip)
        candidate_reason = f"Found candidate D-ring hydroxyl carbon (likely C17) with CIP '{cip}', attached to exocyclic side-chain"
        break  # Stop after the first candidate found

    if candidate_found is not None:
        return True, candidate_reason
    else:
        return False, f"Steroid nucleus identified but no suitable 17α–OH candidate in the D-ring was detected. {candidate_reason}"

# Example usage (for testing – remove or comment when using as module):
if __name__ == "__main__":
    # Test with one example SMILES (17α-hydroxypregnenolone)
    test_smiles = "CC(=O)[C@@]1(O)CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_17alpha_hydroxy_steroid(test_smiles)
    print(result, reason)