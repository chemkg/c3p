"""
Classifies: CHEBI:26255 prenylquinone
"""
#!/usr/bin/env python
"""
Classifies: Prenylquinone

Definition:
    A prenylquinone is defined as a quinone moiety (often a benzoquinone or 1,4–naphthoquinone core)
    substituted by a polyprenyl‐derived side‐chain. Usually, the quinone ring carries two keto groups,
    and the prenyl side‐chain is recognized by an isoprene motif (for example, a “C=C(C)” or “C/C=C(/C)” fragment)
    that is attached outside the quinone core.

This improved classifier attempts:
  1. To detect the quinone motif by first trying two established SMARTS patterns.
  2. If those fail, it scans aromatic six‐membered rings for at least two carbonyl “(C=O)” groups.
  3. It then looks for prenyl fragments using two SMARTS patterns and requires that at least one match is
     found with none of its atoms belonging to the detected quinone core.
"""

from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    The molecule is considered to be a prenylquinone if it contains a quinone core and 
    at least one prenyl (isoprene‐derived) side-chain. For the quinone core, we first try 
    two strict SMARTS patterns (p–benzoquinone and 1,4–naphthoquinone) and if neither is found, 
    we scan aromatic six‐membered rings to see if at least two carbon atoms in the ring have a double‐bonded oxygen.
    For the prenyl side-chain we use two SMARTS patterns: "C=C(C)" and "C/C=C(/C)".
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a prenylquinone, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    quinone_atoms = set()
    found_quinone = False

    # Try two known quinone SMARTS patterns:
    quinone_smarts_list = [
        "O=C1C=CC(=O)C=C1",          # classical p-benzoquinone pattern
        "O=C1C=CC2=C(C=CC(=O)C2=1)"   # 1,4–naphthoquinone core pattern
    ]
    for q_smarts in quinone_smarts_list:
        q_pat = Chem.MolFromSmarts(q_smarts)
        match = mol.GetSubstructMatch(q_pat)
        if match:
            found_quinone = True
            quinone_atoms = set(match)
            break

    # If the SMARTS patterns did not match, try a more generic detection:
    if not found_quinone:
        ring_info = mol.GetRingInfo()
        # Examine each ring (we focus on six-membered aromatic rings, common for quinones)
        for ring in ring_info.AtomRings():
            if len(ring) != 6:
                continue
            # Check if all atoms in the ring are aromatic
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            carbonyl_count = 0
            carbonyl_atom_idxs = set()
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # look for carbon atoms that are double-bonded to an oxygen
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # Check if bond order is double and neighboring atom is oxygen
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                carbonyl_count += 1
                                carbonyl_atom_idxs.add(idx)
                                break
            # If at least 2 carbonyl groups are found in one aromatic ring, assume it is a quinone core.
            if carbonyl_count >= 2:
                found_quinone = True
                quinone_atoms = set(ring)
                break

    if not found_quinone:
        return False, "No quinone motif found"
        
    # Now search for prenyl-derived side-chain(s)
    prenyl_smarts_list = [
        "C=C(C)",     # common isoprene motif
        "C/C=C(/C)"   # alternative stereochemical variant
    ]
    external_prenyl_found = False
    for p_smarts in prenyl_smarts_list:
        p_pat = Chem.MolFromSmarts(p_smarts)
        prenyl_matches = mol.GetSubstructMatches(p_pat)
        for match in prenyl_matches:
            # Check if the prenyl fragment is completely outside the quinone core
            if quinone_atoms.isdisjoint(match):
                external_prenyl_found = True
                break
        if external_prenyl_found:
            break

    if not external_prenyl_found:
        return False, "Quinone core found, but no prenyl side-chain detected outside the quinone ring"

    return True, "Contains a quinone core with a prenyl-derived side-chain"


# For simple testing purposes, you may run:
if __name__ == "__main__":
    # Example test: Napyradiomycin CNQ525.538
    test_smiles = "Br[C@H]1C(O[C@]2(C(=O)C=3C=C(O)C(=C(C3C([C@]2(C1)Cl)=O)O)C)C/C=C(/CCC=C(C)C)\\C)(C)C"
    result, reason = is_prenylquinone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)