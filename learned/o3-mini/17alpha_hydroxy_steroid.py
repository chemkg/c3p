"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI:17alpha-hydroxy steroid
Definition: The alpha-stereoisomer of 17-hydroxy steroid.
Heuristic (improved):
  1. Identify five-membered rings that are fused with at least one six-membered ring.
  2. Build a steroid core as the union of the five-membered ring(s) (candidate “D-ring”) plus
     any neighboring six-membered rings.
  3. For each carbon atom in the steroid core, check if it is chiral and has a direct –OH neighbor
     (an oxygen bearing at least one hydrogen). Also, allow a heavy-atom count (neighbors with atomic
     number > 1) of either 3 or 4. At least one heavy neighbor (ignoring the –OH oxygen) must lie
     outside the steroid core (the exocyclic substituent).
  4. If exactly one candidate is found the molecule is classified as a 17α-hydroxy steroid.
Note: This heuristic is still simplified – stereochemistry and steroid numbering can be complex.
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
    
    # Get all ring information from the molecule.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings detected – not a steroid"
    
    # Collect all six-membered rings (candidate A–C rings) as sets for easy comparison.
    six_membered = [set(ring) for ring in rings if len(ring) == 6]
    
    # Identify candidate D-ring: five-membered rings that are fused (share at least 1 atom) with a six-membered ring.
    candidate_d_rings = []
    for ring in rings:
        if len(ring) == 5:
            ring_set = set(ring)
            # Check if fused with any six-membered ring:
            if any(len(ring_set.intersection(six)) > 0 for six in six_membered):
                candidate_d_rings.append(ring_set)
    
    if not candidate_d_rings:
        return False, "No candidate five-membered D-ring fused with six-membered rings was found."
    
    # Build the fused steroid core: union of the candidate D-rings and any six-membered rings fused to them.
    steroid_core = set()
    for d_ring in candidate_d_rings:
        steroid_core.update(d_ring)
        for six in six_membered:
            if d_ring.intersection(six):
                steroid_core.update(six)
    
    # Optionally sanity-check the size of the steroid core
    if len(steroid_core) < 8:
        return False, "Fused steroid core is too small to be a steroid."
    
    # Now search for the candidate 17α–OH carbon:
    # We look over carbons in the steroid core that are chiral and have an –OH substituent.
    candidates = []
    for atom in mol.GetAtoms():
        # Require carbon atom.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Only consider atoms that are part of the steroid core.
        if atom.GetIdx() not in steroid_core:
            continue
        
        # Check that the carbon is stereochemically defined.
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        
        # Look for an -OH group: an oxygen neighbor that itself has at least one hydrogen.
        oh_found = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                oh_found = True
                break
        if not oh_found:
            continue
        
        # Count heavy atoms (atomic number > 1) among neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # Typical C17 is bound to an -OH, plus two or three carbon neighbors.
        # Allow either 3 or 4 heavy neighbors.
        if len(heavy_neighbors) not in (3, 4):
            continue
        
        # Check that at least one heavy neighbor (excluding the OH oxygen) is exocyclic (i.e. outside the steroid core)
        exocyclic = False
        for nbr in heavy_neighbors:
            if nbr.GetAtomicNum() == 8:
                continue  # skip the OH oxygen
            if nbr.GetIdx() not in steroid_core:
                exocyclic = True
                break
        if not exocyclic:
            continue
        
        # If passed all criteria, record the candidate along with its CIP code if available.
        cip = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else "unknown"
        candidates.append((atom.GetIdx(), cip))
    
    # Evaluate candidate count.
    if len(candidates) == 1:
        return True, f"Found candidate D-ring hydroxyl carbon (likely C17) with CIP '{candidates[0][1]}', attached to exocyclic side-chain"
    elif len(candidates) > 1:
        return False, f"Multiple candidate 17α–OH sites found: {candidates}"
    else:
        return False, "Steroid nucleus identified but no suitable 17α–OH candidate in the D-ring was detected"

# Example usage (for testing – remove or comment when using as a module):
if __name__ == "__main__":
    # Example test: 17α-hydroxypregnenolone
    test_smiles = "CC(=O)[C@@]1(O)CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_17alpha_hydroxy_steroid(test_smiles)
    print(result, reason)