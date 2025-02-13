"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
Definition: Any D-galactoside having beta-configuration at its anomeric centre.
Improvement:
  - We use two explicitly mapped SMARTS patterns that define a beta‐D‐galactopyranoside with
    fixed stereochemistry on the anomeric center and the sequential atoms.
  - Each pattern marks:
      * [*:1]  : the aglycone part (the residue attached via the glycosidic bond)
      * [O:2]  : the glycosidic oxygen
      * [C@H or C@@H:3] : the anomeric carbon (with beta configuration as encoded by the chiral tag)
      * [O:4], [C:5], [C:6], [C:7], [C:8], [C:9] : the remaining atoms in the pyranose ring.
  - Then we check that the ring defined by the mapped ring atoms is a pyranose (6 atoms, exactly 5 carbons and 1 oxygen)
    and that the aglycone atom (mapped as group 1) is not part of that ring.
"""

from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    The method uses two explicitly mapped SMARTS patterns to capture a beta-D-galactopyranoside 
    moiety. After a substructure match, it checks that:
       (i) the matched ring is a pyranose (6-membered ring with five carbons and one oxygen),
      (ii) the anomeric oxygen (of the glycosidic bond) connects the sugar to a residue that is not part 
           of the same pyranose ring (avoiding sugar–sugar linkages).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is identified as a beta-D-galactoside, else False.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define two SMARTS patterns that capture beta-D-galactosides.
    # Pattern 1: Direct glycosidic bond.
    #   [*:1]-[O:2][C@@H:3]1[O:4][C@H:5](O)[C@@H:6](O)[C@@H:7](O)[C@H:8](O)[C@H:9]1O
    # Pattern 2: With a CH2 spacer between the aglycone and the anomeric carbon.
    #   [*:1]-[O:2]C[C@H:3]1[O:4][C@@H:5](O)[C@@H:6](O)[C@H:7](O)[C@H:8](O)[C@H:9]1O
    smarts1 = "[*:1]-[O:2][C@@H:3]1[O:4][C@H:5](O)[C@@H:6](O)[C@@H:7](O)[C@H:8](O)[C@H:9]1O"
    smarts2 = "[*:1]-[O:2]C[C@H:3]1[O:4][C@@H:5](O)[C@@H:6](O)[C@H:7](O)[C@H:8](O)[C@H:9]1O"
    
    patt1 = Chem.MolFromSmarts(smarts1)
    patt2 = Chem.MolFromSmarts(smarts2)
    if patt1 is None or patt2 is None:
        return None, None  # This should not happen if the SMARTS strings are valid

    # Get candidate substructure matches from both patterns
    candidate_matches = mol.GetSubstructMatches(patt1) + mol.GetSubstructMatches(patt2)
    if not candidate_matches:
        return False, "No beta-D-galactoside moiety detected"

    # Get all ring information once
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # List of tuples; each tuple is a set of atom indices in a ring

    # Process each candidate match
    # For our mapped SMARTS, the atom mapping is as follows:
    #   Group 1: aglycone atom (index 0 in match tuple)
    #   Group 2: glycosidic oxygen (index 1)
    #   Group 3: anomeric carbon (index 2)
    #   Groups 4..9: the rest of the sugar ring (indices 3 through 8 in the match tuple)
    for match in candidate_matches:
        # Get the mapped indices from the match tuple
        aglycone_idx = match[0]
        glyco_O_idx   = match[1]
        anomeric_c_idx = match[2]
        
        # The sugar ring is defined by the atoms mapped as 3..9
        sugar_ring_indices = set(match[2:])  # atoms indices corresponding to the ring parts
        
        # Verify that a ring exists which exactly matches the sugar_ring_indices and is 6-membered.
        valid_ring = False
        for ring in atom_rings:
            if set(ring) == sugar_ring_indices and len(ring) == 6:
                # Count atoms: expect exactly 5 carbons and 1 oxygen
                n_C = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                n_O = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                if n_C == 5 and n_O == 1:
                    valid_ring = True
                    break
        if not valid_ring:
            # This candidate does not contain a proper pyranose ring.
            continue

        # Verify that the glycosidic oxygen connects to an atom outside the sugar ring.
        glyco_atom = mol.GetAtomWithIdx(glyco_O_idx)
        external_connection = False
        for nbr in glyco_atom.GetNeighbors():
            # If a neighbor is the aglycone atom and is not part of the sugar ring, consider it external.
            if nbr.GetIdx() == aglycone_idx and nbr.GetIdx() not in sugar_ring_indices:
                external_connection = True
                break
        if not external_connection:
            # In a sugar-sugar linkage the aglycone part would be a sugar atom.
            continue

        # If we have passed all checks, we consider this a valid beta-D-galactoside.
        return True, ("Detected beta-D-galactoside moiety: "
                      "pyranose ring with proper beta stereochemistry and an external glycosidic linkage.")

    # If none of the candidates satisfied all the conditions, return false.
    return False, "Found sugar-like substructure(es), but none met all criteria for a beta-D-galactoside"

# Example test cases
if __name__ == "__main__":
    test_smiles = [
        # cyanidin 3-O-beta-D-galactoside (should pass)
        "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@H]1O",
        # 6-bromo-2-naphthyl beta-D-galactoside (should pass)
        "[C@@H]1([C@@H]([C@@H](O)[C@H]([C@H](O1)CO)O)O)OC=2C=CC=3C=C(C=CC3C2)Br",
        # Methyl beta-D-galactoside (should pass)
        "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",
        # A non-galactoside example (should fail)
        "CC(C)OC(=O)COC1C(O)C(O)C(O)C(O)C1O"
    ]
    for sm in test_smiles:
        flag, reason = is_beta_D_galactoside(sm)
        print(f"SMILES: {sm}\nMatch: {flag}\nReason: {reason}\n")