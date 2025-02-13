"""
Classifies: CHEBI:71971 neoflavonoid
"""
#!/usr/bin/env python
"""
Classifies: neoflavonoid
Definition: Any 1-benzopyran with an aryl substituent at position 4.
Revised heuristic:
  1. Identify candidate pyran rings: six‐membered rings containing exactly one oxygen.
  2. Identify candidate benzene rings: aromatic six‐membered rings.
  3. Find a pyran ring fused (sharing at least two atoms) with a benzene ring.
  4. In that fused pyran ring, examine the atoms that are not in the fused benzene component.
     For each such candidate atom, count external substituents that are part of an independent aromatic ring 
     (one that does not share any atoms with the fused core). If exactly one candidate is found,
     classify the molecule as a neoflavonoid.
If any step fails or we find an ambiguous result (e.g. multiple substituents), the classification is False.
If the SMILES cannot be parsed, return False with an appropriate message.
"""

from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 4.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a neoflavonoid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize and compute aromaticity
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {str(e)}"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # (1) Find candidate pyran rings: 6-membered rings containing exactly one oxygen.
    pyran_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            # Here we require exactly one oxygen in the ring
            if oxy_count == 1:
                pyran_rings.append(set(ring))
                
    if not pyran_rings:
        return False, "No six-membered oxygen-containing (pyran) ring found."
    
    # (2) Find candidate benzene rings: aromatic 6-membered rings.
    benzene_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                benzene_rings.append(set(ring))
    
    if not benzene_rings:
        return False, "No aromatic benzene ring found for fusion."
    
    # (3) Look for a pyran ring fused with a benzene ring (share at least two atoms)
    fused_found = False
    fused_pyran = None
    fused_benzene = None
    for pyran in pyran_rings:
        for benzene in benzene_rings:
            shared = pyran.intersection(benzene)
            if len(shared) >= 2:
                fused_pyran = pyran  # candidate 1-benzopyran core (pyran part)
                fused_benzene = benzene  # fused benzene ring
                fused_found = True
                break
        if fused_found:
            break
    if not fused_found:
        return False, "No fused benzopyran (pyran fused with benzene) core found."
    
    # (4) In a 1-benzopyran, the numbering is such that the position for the aryl substituent (position 4)
    # is on an atom of the pyran that is NOT part of the fused benzene ring.
    # For each such candidate atom, examine neighbors (substituents) outside the entire fused core.
    # We require that exactly one such substituent is found and that it is part of an independent aromatic ring.
    
    candidate_positions = []
    for idx in fused_pyran - fused_benzene:
        atom = mol.GetAtomWithIdx(idx)
        external_aryl_count = 0
        # Examine neighbors that are not in the fused core.
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in fused_pyran:
                continue
            # Check if neighbor is marked aromatic. (This is a preliminary check.)
            if nbr.GetIsAromatic():
                # Now verify if this neighbor is part of an independent aromatic ring.
                # We require that at least one aromatic ring (with 6 or more atoms)
                # contains the neighbor, but none of the ring atoms are in the fused core.
                for ring in atom_rings:
                    if nbr.GetIdx() in ring and len(ring) >= 6:
                        if fused_pyran.intersection(ring) == set():
                            external_aryl_count += 1
                            break
        # Expect exactly one external aryl substituent at the candidate position.
        if external_aryl_count == 1:
            candidate_positions.append(idx)
    
    if len(candidate_positions) == 0:
        return False, "No aryl substituent detected at a candidate position (likely position 4)."
    elif len(candidate_positions) > 1:
        return False, "Multiple candidate positions with aryl substituents detected; ambiguous neoflavonoid assignment."
    
    return True, "Molecule contains a 1-benzopyran core with an aryl substituent at position 4; pattern consistent with a neoflavonoid."

# Example usage:
if __name__ == "__main__":
    # Test with a known prototype: 4-phenylcoumarin (neoflavonoid prototype)
    test_smiles = "C12=CC=CC=C1C(=CC(O2)=O)C3=CC=CC=C3"
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)