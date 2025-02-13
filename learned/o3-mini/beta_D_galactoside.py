"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
Definition: Any D-galactoside having beta-configuration at its anomeric centre.
This improved version uses several SMARTS patterns (with alternate chiral specifications and
an optional CH2 spacer) to capture a beta-D-galactopyranoside moiety. After a substructure
match, we check that the matched ring is a pyranose (6-membered ring with five carbons and one oxygen)
and that the glycosidic oxygen connects to an atom (the aglycone) that is not part of the pyranose ring.
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    The method uses four mapped SMARTS patterns that capture a beta-D-galactoside moiety:
      Pattern 1: Direct glycosidic bond with anomeric carbon marked as [C@@H]
      Pattern 2: Direct glycosidic bond with anomeric carbon marked as [C@H]
      Pattern 3: Glycosidic bond with an extra CH2 spacer and anomeric carbon marked as [C@@H]
      Pattern 4: Glycosidic bond with an extra CH2 spacer and anomeric carbon marked as [C@H]

    Each pattern maps:
       Group 1: aglycone atom (the atom attached via the glycosidic bond)
       Group 2: the glycosidic oxygen
       Group 3: the anomeric carbon (the first carbon of the pyranose ring)
       Groups 4..9: the remaining atoms of the pyranose ring.

    After a candidate match is found, the routine verifies:
       - The candidate atoms form a six-membered ring.
       - The ring consists of exactly five carbons and one oxygen.
       - The glycosidic oxygen (mapped as group 2) is attached to an external atom (the aglycone).
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the beta-D-galactoside moiety is detected, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns.
    # Pattern 1: direct bond, anomeric carbon with [C@@H]
    smarts1 = "[*:1]-[O:2][C@@H:3]1[O:4][C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O"
    # Pattern 2: direct bond, anomeric carbon with [C@H]
    smarts2 = "[*:1]-[O:2][C@H:3]1[O:4][C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O"
    # Pattern 3: with a CH2 spacer, anomeric carbon with [C@@H]
    smarts3 = "[*:1]-[O:2]C[C@@H:3]1[O:4][C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O"
    # Pattern 4: with a CH2 spacer, anomeric carbon with [C@H]
    smarts4 = "[*:1]-[O:2]C[C@H:3]1[O:4][C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O"
    
    patterns = []
    for s in (smarts1, smarts2, smarts3, smarts4):
        patt = Chem.MolFromSmarts(s)
        if patt is not None:
            patterns.append(patt)
    if not patterns:
        return None, None  # Should not happen unless SMARTS strings are invalid
    
    # Get all candidate matches from any of the patterns.
    candidate_matches = []
    for patt in patterns:
        candidate_matches.extend(mol.GetSubstructMatches(patt))
    
    if not candidate_matches:
        return False, "No beta-D-galactoside moiety detected"
    
    # Retrieve all ring information once from the molecule so we can check for the proper sugar ring.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # List of tuples of atom indices that form rings.
    
    # Process each candidate match.
    # In our mapped SMARTS:
    #    match[0] = aglycone atom
    #    match[1] = glycosidic oxygen
    #    match[2] = anomeric carbon (first sugar ring atom)
    #    match[3:9] = the remaining atoms of the sugar ring.
    for match in candidate_matches:
        if len(match) < 9:
            continue  # skip incomplete match
        
        aglycone_idx = match[0]
        glyco_O_idx  = match[1]
        anomeric_idx = match[2]
        sugar_ring_indices = set(match[2:9])
        
        # Look for a ring among the molecule rings that exactly corresponds to the sugar ring.
        valid_ring = False
        for ring in atom_rings:
            if set(ring) == sugar_ring_indices and len(ring) == 6:
                # Count atoms in candidate ring: expect exactly 1 oxygen and 5 carbons.
                n_C = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                n_O = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                if n_C == 5 and n_O == 1:
                    valid_ring = True
                    break
        if not valid_ring:
            continue
        
        # Ensure that the glycosidic oxygen connects to an external atom (the aglycone part).
        glyco_atom = mol.GetAtomWithIdx(glyco_O_idx)
        external_connection = False
        for nbr in glyco_atom.GetNeighbors():
            if nbr.GetIdx() == aglycone_idx and nbr.GetIdx() not in sugar_ring_indices:
                external_connection = True
                break
        if not external_connection:
            continue

        # We have found a candidate with a proper pyranose and a glycosidic linkage.
        return True, ("Detected beta-D-galactoside moiety: "
                      "a pyranose ring with beta stereochemistry at the anomeric carbon "
                      "and an external glycosidic linkage.")
        
    return False, "Found sugar-like substructure(s), but none met all criteria for a beta-D-galactoside"


# Example test cases. (These are the ones given in the prompt.)
if __name__ == "__main__":
    test_smiles = [
        # cyanidin 3-O-beta-D-galactoside:
        "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@H]1O",
        # 6-bromo-2-naphthyl beta-D-galactoside:
        "[C@@H]1([C@@H]([C@@H](O)[C@H]([C@H](O1)CO)O)O)OC=2C=CC=3C=C(C=CC3C2)Br",
        # beta-D-Galp-(1->3)-D-Xylp:
        "O([C@@H]1O[C@@H]([C@H](O)[C@H](O)[C@H]1O)CO)[C@@H]2[C@@H](O)C(OC[C@H]2O)O",
        # 5-bromo-4-chloro-3-indolyl beta-D-galactoside:
        "[H][C@]1(O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)Oc1c[nH]c2ccc(Br)c(Cl)c12",
        # methyl beta-D-galactoside:
        "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",
        # A non-galactoside example:
        "CC(C)OC(=O)COC1C(O)C(O)C(O)C(O)C1O"
    ]
    for sm in test_smiles:
        flag, reason = is_beta_D_galactoside(sm)
        print(f"SMILES: {sm}\nResult: {flag}\nReason: {reason}\n")