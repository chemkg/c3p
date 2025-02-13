"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: 3-oxo-Delta(1) steroid
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
We use heuristic substructure analysis with RDKit to search for:
  1. A ketone group (C(=O)) located in a ring.
  2. Evidence that the ring system is fused (an atom in the ring belongs to at least two rings).
  3. At least one adjacent (alpha) carbon to the ketone that is involved in a non–carbonyl C=C double bond,
     which is our proxy for a Delta(1) double bond in ring A.
Note: This is a heuristic approach and may not capture all edge cases.
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    Criteria:
      - Contains a ketone (C=O) group in a ring (the 3-oxo group),
      - The ring system is fused (at least one atom in the ring belongs to more than one ring),
      - One of the alpha carbons (adjacent to the ketone carbon) is involved in a non-carbonyl C=C double bond,
        corresponding to a double bond between positions 1 and 2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We'll loop over all atoms looking for a candidate ketone carbon.
    # A candidate ketone carbon: atomic number 6; has a double bond to an oxygen atom.
    ring_info = mol.GetRingInfo()
    candidate_found = False
    reason_details = ""
    for atom in mol.GetAtoms():
        # Check for carbon atoms only.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Look for a double bond to oxygen from this carbon (ketone condition)
        ketone = False
        oxy_neighbor = None
        for bond in atom.GetBonds():
            # Check if bond is a double bond
            if bond.GetBondTypeAsDouble() == 2.0:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    ketone = True
                    oxy_neighbor = nbr
                    break
        if not ketone:
            continue
        
        # The candidate ketone carbon must be in a ring (i.e. the "3-oxo" is not exocyclic)
        if not atom.IsInRing():
            # Not in a ring; skip
            continue
        
        # Check that the ring (or one of the rings including this atom) is fused.
        # We determine fusion by checking if the atom belongs to more than one ring
        num_rings_at_atom = ring_info.NumAtomRings(atom.GetIdx())
        fused = num_rings_at_atom > 1
        # Alternatively, if current atom is not fused, check its ring-neighbors.
        if not fused:
            for nbr in atom.GetNeighbors():
                if ring_info.NumAtomRings(nbr.GetIdx()) > 1:
                    fused = True
                    break
        if not fused:
            # Not fused -> likely not part of classical steroid nucleus.
            continue

        # Now check for the C=C (non carbonyl) double bond in one of the alpha carbons.
        # Look at neighbors of the candidate ketone carbon (skip the oxygen we already used)
        delta1_found = False
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == oxy_neighbor.GetIdx():
                continue
            # We expect an alpha carbon (sp2 involvement) within a ring.
            # Check that neighbor is carbon and in a ring.
            if nbr.GetAtomicNum() != 6 or not nbr.IsInRing():
                continue
            
            # Now check if this neighbor participates in a C=C double bond (excluding bonds with the candidate ketone)
            for bond in nbr.GetBonds():
                # Skip the bond connecting nbr and our candidate atom
                if bond.GetOtherAtom(nbr).GetIdx() == atom.GetIdx():
                    continue
                # Check if bond is double and the other atom is carbon.
                if bond.GetBondTypeAsDouble() == 2.0:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 6:
                        delta1_found = True
                        break
            if delta1_found:
                break

        if not delta1_found:
            reason_details = "Candidate 3-oxo group found but no adjacent double bond (Delta(1) candidate) detected."
            continue

        # If we found a fused ketone in a ring with an adjacent double bond, we classify as 3-oxo-Delta(1) steroid.
        candidate_found = True
        break

    if candidate_found:
        return True, "Molecule contains a fused ring ketone (3-oxo group) with an adjacent double bond (Delta(1))."
    else:
        if not reason_details:
            reason_details = "No appropriate 3-oxo group in a fused ring system with an adjacent double bond (Delta(1)) was identified."
        return False, reason_details