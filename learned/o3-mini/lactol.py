"""
Classifies: CHEBI:38131 lactol
"""
#!/usr/bin/env python
"""
Classifies: lactol

Definition:
  Lactols (cyclic hemiacetals) are formed by the intramolecular addition of a hydroxy group to a 
  carbonyl resulting in a cyclic structure where the (former carbonyl) carbon bears two oxygen substituents:
  one coming from an exocyclic hydroxyl (-OH) and one as part of the ring (an acetal linkage). In our heuristic,
  the candidate lactol carbon is required to be in a ring and to be attached to exactly two oxygen atoms:
  one of which is exocyclic and displays at least one explicit hydrogen, and one that is a ring oxygen.
  Furthermore, the lactol carbon and its ring oxygen neighbor must belong together to a 5- or 6-membered ring.
  
Note:
  This improved heuristic attempts to reduce false positives (e.g. in large glycosylated structures)
  while still detecting common lactol motifs (including simple sugars like alpha-D-fructopyranose)
  and more complex natural products. It may still fail on unusual structures.
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol (cyclic hemiacetal) based on its SMILES string.
    
    A lactol (cyclic hemiacetal) is produced by the intramolecular addition of a hydroxy group to a carbonyl 
    (aldehydic or ketonic) center, yielding a ring in which that carbon (the candidate lactol carbon) is bound 
    to exactly two oxygen atoms: one as part of the ring (forming an acetal linkage) and one as an exocyclic -OH.
    In this heuristic:
      - The candidate lactol carbon (a sp3 carbon in a ring) must have exactly two heavy-atom oxygen neighbors.
      - Out of these, one oxygen must be in a ring and the other must be exocyclic and show a directly attached hydrogen.
      - In addition, the candidate carbon and its ring oxygen neighbor must share a ring of size 5 or 6.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lactol, False otherwise.
        str: Reason for the decision.
    """
    # Parse SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (so that -OH hydrogens are explicit)
    mol = Chem.AddHs(mol)
    
    # Get ring information (a list of atom-index tuples for each ring)
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Loop through all atoms to search for a candidate lactol carbon.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # must be carbon
            continue
        if not atom.IsInRing():
            continue
        
        # Get all oxygen neighbors (heavy atom oxygen)
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        # We now require that the lactol carbon be attached to exactly two oxygens.
        if len(oxy_neighbors) != 2:
            continue
        
        # Separate the oxygen neighbors into (a) ring-bound and (b) exocyclic -OH candidates.
        ring_oxy_list = []
        exo_oh_list = []
        for oxy in oxy_neighbors:
            if oxy.IsInRing():
                ring_oxy_list.append(oxy)
            else:
                # Check that this exocyclic oxygen bears at least one hydrogen (–OH)
                if any(n.GetAtomicNum() == 1 for n in oxy.GetNeighbors()):
                    exo_oh_list.append(oxy)
        
        # We require exactly one ring oxygen and exactly one exocyclic -OH oxygen.
        if len(ring_oxy_list) != 1 or len(exo_oh_list) != 1:
            continue
        
        # Now, check that the candidate carbon and the ring oxygen neighbor are together part of some 5- or 6-membered ring.
        candidate_idx = atom.GetIdx()
        ring_oxy_idx = ring_oxy_list[0].GetIdx()
        ring_match_found = False
        for ring in ring_info:
            if candidate_idx in ring and ring_oxy_idx in ring:
                if len(ring) in (5, 6):
                    ring_match_found = True
                    break
        
        if ring_match_found:
            return True, "Found lactol moiety: candidate carbon in a 5- or 6-membered ring with one exocyclic OH and one ring oxygen"
    
    # If no candidate meets the criteria
    return False, "No suitable lactol moiety pattern detected in the molecule"

# Example usage (can be removed when integrated into a larger pipeline)
if __name__ == "__main__":
    test_smiles = [
        "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha-D-fructopyranose (should be classified as lactol)
        "O=C1OC(CO)C(O)C1"  # a cyclic structure that may not fulfill the strict criteria
    ]
    for s in test_smiles:
        result, reason = is_lactol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")