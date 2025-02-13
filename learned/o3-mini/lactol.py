"""
Classifies: CHEBI:38131 lactol
"""
#!/usr/bin/env python
"""
Classifies: lactol
Definition:
  Lactols (cyclic hemiacetals) are formed by the intramolecular addition of a hydroxy group
  to an aldehydic or ketonic carbonyl group. In the resulting structure (1-oxacycloalkan-2-ol or 
  unsaturated analogues) the (former) carbonyl carbon bears two oxygen substituents: one belonging 
  to the ring (forming an acetal linkage) and one as an exocyclic hydroxyl.
  
This program uses a heuristic approach:
  - We add explicit hydrogens so that hydroxyl groups are properly recognized.
  - We loop over carbon atoms that are in rings.
  - For each such carbon, we check whether it has at least one exocyclic hydroxyl group 
    (an oxygen not in a ring that has a hydrogen) and at least one oxygen neighbor that is in a ring.
  - We then verify that the carbon and its ring oxygen neighbor are part of a 5- or 6-membered ring.
  
Note: this heuristic may still miss some lactol cases or may misclassify very unusual structures.
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol (cyclic hemiacetal) based on its SMILES string.
    
    A lactol results from the intramolecular addition of a hydroxy group to a carbonyl,
    yielding a cyclic structure where the (former carbonyl) carbon bears both an exocyclic
    hydroxyl group and an oxygen that is part of the ring.
    
    Heuristic criteria:
      - The candidate lactol carbon must be in a ring.
      - It must have at least one exocyclic oxygen atom that has at least one hydrogen 
        (i.e. -OH) and one oxygen neighbor that is part of a ring.
      - The carbon and its ring oxygen neighbor must belong to a 5- or 6-membered ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lactol, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydroxyl hydrogens are represented
    mol = Chem.AddHs(mol)
    
    # Get ring information (a list of tuples, each is a ring defined by atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Loop over all atoms. We are interested in candidate lactol carbons.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # candidate carbon
            continue
        if not atom.IsInRing():
            continue

        # Get all neighbors that are oxygen atoms.
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(oxy_neighbors) < 2:
            # Lactol carbon must be bound to at least 2 oxygens (one from -OH and one ring oxygen)
            continue

        exo_oh_found = False
        ring_oxy_found = None

        # Check each oxygen neighbor to see if one is an exocyclic hydroxyl (-OH not in a ring)
        # and one is an oxygen that is in a ring.
        for oxy in oxy_neighbors:
            if oxy.IsInRing():
                # save one candidate that is in the ring
                ring_oxy_found = oxy
            else:
                # exocyclic candidate: it should have at least one hydrogen explicit
                has_h = any(nbr.GetAtomicNum() == 1 for nbr in oxy.GetNeighbors())
                if has_h:
                    exo_oh_found = True

        if not (exo_oh_found and ring_oxy_found):
            # We need both an exocyclic -OH and a ring-bound oxygen.
            continue
        
        # Now ensure that the candidate carbon and the ring oxygen neighbor are part of a 
        # ring of size 5 or 6. Check each ring that includes the candidate carbon.
        lactol_ring_found = False
        carbon_idx = atom.GetIdx()
        ring_oxy_idx = ring_oxy_found.GetIdx()
        for ring in ring_info:
            if carbon_idx in ring and ring_oxy_idx in ring:
                if len(ring) in [5, 6]:
                    lactol_ring_found = True
                    break
        
        if lactol_ring_found:
            return True, "Found lactol moiety: candidate carbon in a 5- or 6-membered ring with an exocyclic OH and ring oxygen neighbor"
    
    return False, "No suitable lactol moiety pattern detected in the molecule"

# Example usage (can be removed when integrated in a larger pipeline):
if __name__ == "__main__":
    test_smiles = [
        "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha-D-fructopyranose (valid cyclic hemiacetal)
        "O=C1OC(CO)C(O)C1"  # A simple ring; may or may not be considered lactol
    ]
    for s in test_smiles:
        result, reason = is_lactol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")