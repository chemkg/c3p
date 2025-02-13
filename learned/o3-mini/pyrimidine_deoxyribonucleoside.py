"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: A deoxyribonucleoside containing a pyrimidine base.
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside.
    A pyrimidine deoxyribonucleoside contains a deoxyribose sugar (a furanose ring missing an OH at the 2' position)
    that is linked via a glycosidic bond to a pyrimidine base.
    
    This function uses heuristic substructure searches:
      - It searches for a pyrimidine ring using a generic SMARTS "c1cncnc1".
      - It searches for a candidate deoxyribose sugar by looking for a five-membered ring 
        that has exactly one oxygen in the ring and that has three exocyclic –OH (or CH2OH-like) attachments.
      - It then confirms if one or more atoms in the sugar ring is connected to an atom that is part of 
        the pyrimidine base.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a pyrimidine deoxyribonucleoside, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Search for a pyrimidine base ---
    # Use a generic pyrimidine pattern: a six-membered aromatic ring with alternating C and N.
    pyrimidine_smarts = "c1cncnc1"
    pyrimidine_query = Chem.MolFromSmarts(pyrimidine_smarts)
    if not mol.HasSubstructMatch(pyrimidine_query):
        return False, "Pyrimidine base not found"
    
    # Get atom indices from one pyrimidine match
    pyr_matches = mol.GetSubstructMatches(pyrimidine_query)
    # Take the first match (if multiple found, one candidate is enough)
    pyr_atom_indices = set(pyr_matches[0])
    
    # --- Step 2. Search for a deoxyribose sugar moiety ---
    # Our approach: Look for a five-membered ring containing exactly one oxygen (i.e. a furanose ring).
    # Then, count –OH (or CH2OH) attachments on the ring atoms.
    ring_info = mol.GetRingInfo()
    candidate_sugar = None
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Count number of ring atoms that are oxygen
            oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O")
            if oxy_in_ring == 1:
                # Count exocyclic hydroxyl (or –OH/CH2OH) attachments from ring atoms.
                # We consider an exocyclic oxygen if it is not part of the ring and it has at least one hydrogen.
                exo_oh_count = 0
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() not in ring and nbr.GetSymbol() == "O":
                            # Check if this oxygen has at least one hydrogen attached.
                            if any(neigh.GetAtomicNum() == 1 for neigh in nbr.GetNeighbors()):
                                exo_oh_count += 1
                # For a deoxyribose sugar (in a nucleoside context) we expect three oxygen attachments:
                # one CH2OH (typically the 5'-group) and two -OH groups (typically at the 3' and, if present, at positions not 2').
                if exo_oh_count == 3:
                    # Finally, verify that the sugar ring is connected to the pyrimidine base.
                    connected = False
                    for idx in ring:
                        atom = mol.GetAtomWithIdx(idx)
                        for nbr in atom.GetNeighbors():
                            if nbr.GetIdx() in pyr_atom_indices:
                                connected = True
                                break
                        if connected:
                            break
                    if connected:
                        candidate_sugar = ring
                        break  # Found an appropriate sugar ring

    if candidate_sugar is None:
        return False, "Deoxyribose sugar moiety not found"

    return True, "Molecule contains a pyrimidine base linked to a deoxyribose sugar"

# For testing purposes, you can try calling the function with one of the provided SMILES
if __name__ == "__main__":
    # Example: 2'-deoxyuridine
    test_smiles = "OC[C@H]1O[C@H](C[C@@H]1O)n1ccc(=O)[nH]c1=O"
    result, reason = is_pyrimidine_deoxyribonucleoside(test_smiles)
    print(result, reason)