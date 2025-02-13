"""
Classifies: CHEBI:26255 prenylquinone
"""
#!/usr/bin/env python
"""
Classifies: Prenylquinone

Definition:
    A prenylquinone is a quinone moiety (often a p-benzoquinone, 1,4–naphthoquinone or fused variant)
    substituted by one or more polyprenyl-derived side-chains.
    
    This classifier:
      1. Tries two strict quinone SMARTS patterns (for p-benzoquinone and 1,4–naphthoquinone).
      2. If those fail, scans all aromatic rings (of length 6–10) for at least two carbonyl (C=O) bonds.
      3. Then searches for prenyl fragments using two prenyl SMARTS patterns and requires that at least one
         such fragment is found that is entirely (a) outside the detected quinone core and (b) not part of any ring.
"""

from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    The molecule is considered prenylquinone if:
      - It contains a quinone core. This is established either by a strict match of one of two SMARTS patterns
        or by scanning an aromatic ring (of length 6-10) for at least 2 carbonyl (C=O) groups.
      - It contains at least one prenyl (isoprene-derived) side-chain fragment that is outside both any ring and 
        the quinone core.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple with True and a reason if classified as prenylquinone;
                   else False and a reason describing the shortcoming.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found_quinone = False
    quinone_atoms = set()
    
    # First try two established quinone SMARTS patterns.
    quinone_smarts_list = [
        "O=C1C=CC(=O)C=C1",         # classical p-benzoquinone pattern
        "O=C1C=CC2=C(C=CC(=O)C2=1)"   # 1,4–naphthoquinone core pattern
    ]
    for q_smarts in quinone_smarts_list:
        q_pat = Chem.MolFromSmarts(q_smarts)
        match = mol.GetSubstructMatch(q_pat)
        if match:
            found_quinone = True
            quinone_atoms = set(match)
            break

    # If not found by SMARTS, try generic detection by scanning aromatic rings.
    if not found_quinone:
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if not (6 <= len(ring) <= 10):
                continue
            # Check if all atoms in the ring are aromatic.
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            carbonyl_count = 0
            # Check for at least 2 C=O bonds in the ring.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # only care about carbon atoms
                    for bond in atom.GetBonds():
                        # Check if bond is double and the neighbor is oxygen.
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
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

    # Now search for prenyl-derived side chains.
    # We use two SMARTS patterns that capture isoprene-like fragments.
    prenyl_smarts_list = [
        "C=C(C)",     # common isoprene motif
        "C/C=C(/C)"   # alternative stereochemical variant
    ]
    external_prenyl_found = False
    for p_smarts in prenyl_smarts_list:
        p_pat = Chem.MolFromSmarts(p_smarts)
        prenyl_matches = mol.GetSubstructMatches(p_pat)
        for match in prenyl_matches:
            # Check that the prenyl fragment does not overlap with quinone core.
            if not quinone_atoms.isdisjoint(match):
                continue
            # Also require that none of the atoms in this prenyl match are part of any ring.
            if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
                continue
            external_prenyl_found = True
            break
        if external_prenyl_found:
            break

    if not external_prenyl_found:
        return False, "Quinone core found, but no prenyl side-chain detected outside rings and outside the quinone core"

    return True, "Contains a quinone core with a prenyl-derived side-chain"

# For basic testing purposes:
if __name__ == "__main__":
    # Example test: Napyradiomycin CNQ525.538 (one of the real prenylquinones)
    test_smiles = "Br[C@H]1C(O[C@]2(C(=O)C=3C=C(O)C(=C(C3C([C@]2(C1)Cl)=O)O)C)C/C=C(/CCC=C(C)C)\\C)(C)C"
    result, reason = is_prenylquinone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)