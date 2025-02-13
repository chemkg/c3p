"""
Classifies: CHEBI:26255 prenylquinone
"""
#!/usr/bin/env python
"""
Classifies: Prenylquinone

Definition:
    A prenylquinone is a quinone moiety (often a p-benzoquinone, 1,4–naphthoquinone or fused variant)
    substituted by one or more polyprenyl-derived side-chains.
    
    In this revised classifier:
      1. We search for a quinone core using two SMARTS patterns (p-benzoquinone and an extended 1,4-naphthoquinone)
         and, if necessary, a generic search scanning rings (size 6–10) for at least two carbonyl (C=O) groups.
      2. We then look for prenyl side-chain fragments using two prenyl SMARTS patterns.
         In addition, we require that at least one such prenyl fragment:
           (a) is attached (via a carbon–carbon bond) directly to the quinone core,
           (b) has some atoms not in any ring (i.e. extends outward).
           
Note:
    The algorithm is still heuristic. Improvements included expanding the quinone SMARTS,
    relaxing the aromaticity requirement in generic detection, and adding tests on the prenyl fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    The molecule is considered a prenylquinone if:
      - It contains a quinone core. This is established either by a strict match of one of two SMARTS patterns
        or by scanning rings (of length 6–10) for at least 2 carbonyl (C=O) groups.
      - It contains at least one prenyl (isoprene-derived) side-chain fragment that:
            * is attached via a C–C bond to the quinone core, and
            * extends outside of any ring.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple with True and a reason if classified as prenylquinone;
                   else False and a reason describing the shortcoming.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    found_quinone = False
    quinone_atoms = set()
    
    # Step 1: Look for the quinone core using established SMARTS patterns.
    quinone_smarts_list = [
        "O=C1C=CC(=O)C=C1",           # classical p-benzoquinone pattern
        "O=C1C=CC2=C(C=CC(=O)C2)C1"     # extended 1,4-naphthoquinone fused pattern
    ]
    for q_smarts in quinone_smarts_list:
        q_pat = Chem.MolFromSmarts(q_smarts)
        match = mol.GetSubstructMatch(q_pat)
        if match:
            found_quinone = True
            quinone_atoms = set(match)
            break

    # If strict patterns failed, try generic detection.
    if not found_quinone:
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if not (6 <= len(ring) <= 10):
                continue
            # Count carbonyl groups (C=O double bonds) in the ring.
            carbonyl_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # consider only carbon atoms in the ring
                    for bond in atom.GetBonds():
                        # Check for a double bond to oxygen.
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                carbonyl_count += 1
                                break
            if carbonyl_count >= 2:
                found_quinone = True
                quinone_atoms = set(ring)
                break

    if not found_quinone:
        return False, "No quinone motif found"
    
    # Step 2: Search for prenyl-derived side chains.
    # We use two SMARTS patterns that capture isoprene-like fragments.
    prenyl_smarts_list = [
        "[C;!R]=[C;!R]([CH3])",  # basic isoprene unit (non-ring atoms)
        "C/C=C(/C)C"             # extended stereochemical variant
    ]
    prenyl_found = False
    for p_smarts in prenyl_smarts_list:
        p_pat = Chem.MolFromSmarts(p_smarts)
        prenyl_matches = mol.GetSubstructMatches(p_pat)
        for match in prenyl_matches:
            # Skip very short matches.
            if len(match) < 3:
                continue
            
            # Check that at least one atom in the prenyl match is directly bonded
            # (via a C–C bond) to an atom in the quinone core.
            attached_to_quinone = False
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in quinone_atoms:
                        if atom.GetAtomicNum() == 6 and nbr.GetAtomicNum() == 6:
                            attached_to_quinone = True
                            break
                if attached_to_quinone:
                    break
            if not attached_to_quinone:
                continue
            
            # Require that the prenyl fragment extends outside any ring.
            if all(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
                continue
            
            prenyl_found = True
            break
        if prenyl_found:
            break

    if not prenyl_found:
        return False, "Quinone core found, but no valid prenyl side-chain detected attached to it"
    
    return True, "Contains a quinone core with a prenyl-derived side-chain"

# For basic testing purposes:
if __name__ == "__main__":
    # Example test: Napyradiomycin CNQ525.538 (one of the real prenylquinones)
    test_smiles = "Br[C@H]1C(O[C@]2(C(=O)C=3C=C(O)C(=C(C3C([C@]2(C1)Cl)=O)O)C)C/C=C(/CCC=C(C)C)\\C)(C)C"
    result, reason = is_prenylquinone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)