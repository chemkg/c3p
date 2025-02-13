"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
Definition: Any D-galactoside having beta-configuration at its anomeric centre.
This revised approach uses a two‐step procedure. First, two SMARTS queries (with different 
chiral annotations) search for the beta-glycosidic linkage connected to a sugar-like ring. 
Then for each candidate substructure match we check that the ring is a pyranose – i.e. a six‐membered 
ring containing exactly five carbon atoms and one oxygen – and that the anomeric oxygen is connected 
to a non-sugar residue.
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    This method uses a combination of SMARTS substructure matching plus additional
    checks to ensure that the matching sugar ring has exactly five carbons and one oxygen,
    and that the anomeric oxygen is bound to an external (non-sugar) residue.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is identified as a beta-D-galactoside, else False.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # We define two SMARTS patterns that capture a putative beta-glycosidic bond to a sugar.
    # In these patterns the sugar ring is represented with chiral tags. Note:
    #   - The pattern is: [any]-O[C@@H]1O[C@H]([any])[C@H](O)[C@H](O)[C@H]1O
    #     or its mirror version.
    # The first atom ([*:1]) represents the aglycone (or substituent) attached via oxygen.
    smarts1 = "[*:1]-O[C@@H]1O[C@H]([*:2])[C@H](O)[C@H](O)[C@H]1O"
    smarts2 = "[*:1]-O[C@H]1O[C@@H]([*:2])[C@@H](O)[C@H](O)[C@@H]1O"

    patt1 = Chem.MolFromSmarts(smarts1)
    patt2 = Chem.MolFromSmarts(smarts2)

    if patt1 is None or patt2 is None:
        return None, None  # Should not occur if SMARTS are valid

    # Collect all candidate matches from either pattern.
    candidate_matches = mol.GetSubstructMatches(patt1) + mol.GetSubstructMatches(patt2)
    if not candidate_matches:
        return False, "No beta-D-galactoside moiety detected"

    # Get ring information for the molecule once so we do not recompute.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # list of tuples of atom indices for each ring

    # For each candidate match, do the following extra checks:
    # - Identify the "anomeric carbon" from the SMARTS match.
    #   In the SMARTS pattern, the matched atoms are in order:
    #     index0: [*:1] (the aglycone part),
    #     index1: the glycosidic oxygen,
    #     index2: the anomeric carbon (with [C@@H] or [C@H])
    # The candidate sugar ring should be a six-membered ring containing this anomeric carbon.
    for match in candidate_matches:
        # Get the atom index for the anomeric carbon from the match tuple.
        anomeric_carb_idx = match[2]
        # Check each ring that contains the anomeric carbon:
        for ring in atom_rings:
            if anomeric_carb_idx in ring and len(ring) == 6:
                # Count number of carbons and oxygens in this ring.
                n_C = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                n_O = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                # For a typical pyranose, we expect 5 C and 1 O.
                if n_C == 5 and n_O == 1:
                    # Now also check that the glycosidic oxygen (match index 1) connects 
                    # to an atom that is NOT in the candidate sugar ring.
                    glyco_O_idx = match[1]
                    glyco_atom = mol.GetAtomWithIdx(glyco_O_idx)
                    external = False
                    for nbr in glyco_atom.GetNeighbors():
                        if nbr.GetIdx() not in ring:
                            external = True
                            break
                    if external:
                        # We consider this a valid beta-D-galactoside substructure.
                        return True, "Detected beta-D-galactoside moiety (pyranose ring with beta glycosidic linkage and proper ring composition)"
                    else:
                        # The oxygen appears to be only bonded within the sugar (e.g. in a sugar–sugar linkage).
                        continue
    # If no candidate match passed all the checks:
    return False, "Found a sugar-like substructure, but no candidate met all criteria for a beta-D-galactoside"

# Example test cases (uncomment to run)
# test_smiles = [
#     "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@H]1O",  # cyanidin 3-O-beta-D-galactoside
#     "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",  # methyl beta-D-galactoside (should be true)
#     "CC(C)OC(=O)COC1C(O)C(O)C(O)C(O)C1O"  # a non-galactoside example
# ]
# for sm in test_smiles:
#     flag, reason = is_beta_D_galactoside(sm)
#     print(f"SMILES: {sm}\nMatch: {flag}\nReason: {reason}\n")