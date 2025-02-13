"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: neoflavonoid 
Definition: Any 1-benzopyran with an aryl substituent at position 4.
The heuristic implemented below uses the following steps:
  1. Identify candidate pyran rings: rings with 6 atoms that contain exactly one oxygen.
  2. Identify candidate benzene rings: aromatic 6-membered rings.
  3. Check for a pyran ring fused with a benzene ring (they must share at least 2 atoms).
  4. In the pyran ring, look for an atom (the candidate “position 4”) that is not in the fused benzene ring,
     and that carries a substituent (outside the pyran ring) which is an aromatic ring (of 6 or more atoms).
If such a pattern is found, we return True with an explanation.
If not, we return False with a reason.
If the input SMILES is not valid, we return False with an “Invalid SMILES string” message.
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
        
    # Ensure aromaticity is computed
    Chem.SanitizeMol(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Find candidate pyran rings: 6-membered rings containing exactly one oxygen.
    pyran_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            # Allow exactly one oxygen in the ring
            if oxy_count == 1:
                pyran_rings.append(set(ring))
    
    if not pyran_rings:
        return False, "No six-membered oxygen-containing (pyran) ring found."
    
    # Find candidate benzene rings: aromatic six-membered rings.
    benzene_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            # Check all atoms in ring are aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                benzene_rings.append(set(ring))
    
    if not benzene_rings:
        return False, "No aromatic benzene ring found for fusion."
    
    # Look for a candidate pyran ring fused with a benzene ring (share at least two atoms)
    fused_found = False
    for pyran in pyran_rings:
        for benzene in benzene_rings:
            shared = pyran.intersection(benzene)
            if len(shared) >= 2:
                fused_pyran = pyran  # our candidate benzopyran core
                fused_benzene = benzene
                fused_found = True
                break
        if fused_found:
            break
    if not fused_found:
        return False, "No fused benzopyran (pyran fused with benzene) core found."
    
    # In a 1-benzopyran the numbering is such that one carbon in the pyran (not part of the benzene fusion)
    # should carry an extra substituent that is an aryl (benzene-like) ring.
    # Check each atom in the pyran that is not in the fused benzene ring for an external aryl substituent.
    candidate_position_found = False
    for idx in fused_pyran - fused_benzene:
        atom = mol.GetAtomWithIdx(idx)
        # Look at neighbors that are not in the pyran ring:
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in fused_pyran:
                continue
            # Check if the neighbor is aromatic:
            if nbr.GetIsAromatic():
                # Check if this neighbor is part of an independent aromatic ring of 6 atoms
                for ring in atom_rings:
                    if nbr.GetIdx() in ring and len(ring) >= 6:
                        # Ensure the ring is not simply a small substituent fragment by checking that
                        # not too many of its atoms are in the benzopyran core.
                        if len(set(ring).intersection(fused_pyran)) < 3:
                            candidate_position_found = True
                            break
                if candidate_position_found:
                    break
        if candidate_position_found:
            break
    
    if not candidate_position_found:
        return False, "No aryl substituent detected at a candidate position (likely position 4)."
    
    return True, "Molecule contains a 1-benzopyran core with an aryl substituent at position 4; pattern consistent with a neoflavonoid."

# Below are examples on how to call the function.
if __name__ == "__main__":
    # Example SMILES for 4-phenylcoumarin, a neoflavonoid prototype:
    test_smiles = "C12=CC=CC=C1C(=CC(O2)=O)C3=CC=CC=C3"
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)