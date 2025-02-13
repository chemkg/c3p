"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: Beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
This implementation uses a combined substructure and stereochemical “ring‐walking” approach:
  1. It first assigns stereochemistry.
  2. It then scans for 6-membered rings (pyranoses) that have exactly one oxygen (ring atom)
     and five carbon atoms.
  3. For each such ring the algorithm iterates over its carbon atoms and looks for a candidate
     anomeric carbon – the candidate must have at least one exocyclic oxygen (i.e. not in the ring)
     and its CIP code (as assigned by RDKit) must be “R” (a proxy for the beta configuration).
  4. If such a candidate is found, the ring is “walked” in a defined order: the candidate anomeric
     carbon (C1) should be bonded to two ring neighbours. One of those is expected to carry the –CH2OH
     substituent (i.e. the sugar’s C5). With a simple traversal we then obtain an order for four ring carbons
     (C1–C4) whose CIP codes are compared with the expected (R, S, R, R) for beta-D-glucopyranose.
     
If a ring meets all these criteria the molecule is classified as a beta-D-glucoside.
Due to the diversity of SMILES representations and sugar substitutions this algorithm is heuristic.
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    It does so by scanning for a pyranose (6-membered) ring that meets:
      - The ring has exactly 1 oxygen and 5 carbons.
      - One of the ring carbons (candidate anomeric carbon) is directly bonded to an exocyclic oxygen.
      - That candidate has defined stereochemistry with a CIP code of "R" (as expected for beta anomers).
      - Using the ring connectivity we “order” the five ring carbons (which in a D-glucopyranose
        should be numbered C1–C5, with C1 = anomeric, and C2, C3, C4 having expected CIP values).
        The expected CIP fingerprints (when traversed in order) for beta-D-glucose are:
          C1: "R", C2: "S", C3: "R", C4: "R"
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if a beta-D-glucoside fragment is detected, False otherwise.
        str: Explanation for the classification result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Assign stereochemistry (force and clean)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    mol.UpdatePropertyCache()
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Expected stereochemical fingerprint for beta-D-glucopyranose.
    expected_CIP = ("R", "S", "R", "R")  # for C1, C2, C3, C4 (C5 carries CH2OH)
    
    # Look over each ring that is 6-membered.
    for ring in rings:
        if len(ring) != 6:
            continue
        # Count atoms in ring: we expect 1 oxygen and 5 carbons.
        o_count = 0
        c_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                o_count += 1
            elif atom.GetAtomicNum() == 6:
                c_count += 1
        if o_count != 1 or c_count != 5:
            continue
        
        # Loop over ring atoms to search for a candidate anomeric carbon.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Candidate must be carbon.
            if atom.GetAtomicNum() != 6:
                continue
            # Look for an exocyclic oxygen neighbor (i.e. not in the ring).
            exo_oxys = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring]
            if not exo_oxys:
                continue
            # Check that stereochemistry is specified and CIP code is "R"
            try:
                cip = atom.GetProp("_CIPCode")
            except KeyError:
                # CIP code not set -- skip this candidate.
                continue
            if cip != "R":
                continue
            # At this point, we treat this atom as the candidate anomeric carbon (C1).
            candidate_C1 = atom
            # Its ring neighbors (should be exactly two)
            ring_neigh_idxs = [nbr.GetIdx() for nbr in candidate_C1.GetNeighbors() if nbr.GetIdx() in ring]
            if len(ring_neigh_idxs) != 2:
                continue
            # We expect one of these neighbors to have an exocyclic CH2OH (the sugar C5).
            candidate_C5 = None
            candidate_C2 = None
            for nb_idx in ring_neigh_idxs:
                nb = mol.GetAtomWithIdx(nb_idx)
                # Look for an exocyclic carbon (not in ring) attached to this neighbor having exactly 2 hydrogens.
                for nb2 in nb.GetNeighbors():
                    if nb2.GetIdx() in ring:
                        continue
                    if nb2.GetAtomicNum() == 6 and nb2.GetTotalNumHs() == 2:
                        candidate_C5 = nb
                        break
                if candidate_C5 is not None:
                    break
            # If we did not find such a neighbor, then try the other possibility:
            if candidate_C5 is None:
                continue
            # Then the other ring neighbor of candidate C1 becomes candidate C2.
            for nb_idx in ring_neigh_idxs:
                if nb_idx != candidate_C5.GetIdx():
                    candidate_C2 = mol.GetAtomWithIdx(nb_idx)
                    break
            if candidate_C2 is None:
                continue
            # Now, traverse the ring to obtain an order: expected order: C1 (candidate), then C2,
            # then C3 the unique ring neighbor of C2 (other than C1), then C4 the unique neighbor of C3 (other than C2),
            # and finally expecting candidate_C5 as the other neighbor of C4.
            # Get C3:
            neighbors_C2 = [nbr for nbr in candidate_C2.GetNeighbors() if nbr.GetIdx() in ring and nbr.GetIdx() != candidate_C1.GetIdx() and nbr.GetAtomicNum()==6]
            if len(neighbors_C2) != 1:
                continue
            candidate_C3 = neighbors_C2[0]
            # Get C4 from candidate_C3 (the ring neighbor not equal to candidate_C2).
            neighbors_C3 = [nbr for nbr in candidate_C3.GetNeighbors() if nbr.GetIdx() in ring and nbr.GetIdx() != candidate_C2.GetIdx() and nbr.GetAtomicNum()==6]
            if len(neighbors_C3) != 1:
                continue
            candidate_C4 = neighbors_C3[0]
            # Now, candidate_C4 should be ring-connected to candidate_C5.
            ring_neigh_C4 = [nbr for nbr in candidate_C4.GetNeighbors() if nbr.GetIdx() in ring and nbr.GetAtomicNum()==6]
            if candidate_C5.GetIdx() not in [a.GetIdx() for a in ring_neigh_C4]:
                continue
            
            # Check that all four candidate carbons have defined CIP codes.
            try:
                cip_C1 = candidate_C1.GetProp("_CIPCode")
                cip_C2 = candidate_C2.GetProp("_CIPCode")
                cip_C3 = candidate_C3.GetProp("_CIPCode")
                cip_C4 = candidate_C4.GetProp("_CIPCode")
            except KeyError:
                continue
            # Compare with expected fingerprint (for beta-D-glucose: C1=R, C2=S, C3=R, C4=R).
            observed = (cip_C1, cip_C2, cip_C3, cip_C4)
            if observed == expected_CIP:
                return True, "Beta-D-glucoside fragment detected with correct ring stereochemistry."
    return False, "No beta-D-glucoside fragment found in the molecule."

# Example usage for testing:
if __name__ == "__main__":
    test_examples = [
        # One example that should be classified as beta-D-glucoside:
        ("OC[C@H]1O[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "Example with a beta-D-glucoside fragment"),
        # You can add additional SMILES strings for testing.
    ]
    for smi, desc in test_examples:
        result, reason = is_beta_D_glucoside(smi)
        print("SMILES:", smi)
        print("Description:", desc)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)