"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: CHEBI:3-oxo-Δ(1) steroid
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
Heuristic approach:
  1. Find a candidate six-membered ring carbon bearing a ketone (C(=O)).
  2. Check that this ketone carbon is connected to another six-membered ring carbon that participates in a C=C double bond (as a proxy for a Delta(1) bond).
  3. Ensure that the candidate ketone is in a fused ring system (i.e. at least one atom in the candidate fragment appears in more than one ring) and that the overall molecule has at least 3 rings.
Note: This heuristic does not cover all edge cases.
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Δ(1) steroid based on its SMILES string.
    
    Criteria (heuristic):
      - Has at least 3 rings (typical steroids have 4 fused rings),
      - Contains a six-membered ring carbon with a ketone group (C(=O)),
      - That ketone carbon is directly linked to a six-membered ring carbon which
        participates in a non-carbonyl C=C double bond (proxy for the Δ(1) double bond),
      - And at least one atom of the candidate fragment is in a fused ring (belongs to >1 ring).
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if the molecule qualifies as a 3-oxo-Δ(1) steroid, else False.
       str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    # Steroid nucleus usually contains 3 or 4 rings; we require at least 3 rings.
    if ring_info.NumRings() < 3:
        return False, "Molecule does not have enough rings (<3) to be a steroid"

    # SMARTS pattern explanation:
    #   [$([#6;R6](=O))]   : a six-membered ring carbon (atomic number 6) that has a double bond to oxygen (ketone).
    #   -                   : a single bond to
    #   [$([#6;R6])]       : another six-membered ring carbon,
    #   =                   : which is double-bonded to
    #   [$([#6;R6])]       : another six-membered ring carbon.
    pattern_smarts = "[$([#6;R6](=O))]-[$([#6;R6])]=[$([#6;R6])]"
    patt = Chem.MolFromSmarts(pattern_smarts)
    if patt is None:
        return False, "Error in SMARTS pattern"
    
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return False, "No six-membered ring ketone adjacent to a ring double bond (Δ(1) candidate) found"

    # For each candidate match, check if at least one of the atoms in the match is fused (i.e. belongs to >1 ring)
    for match in matches:
        fused = False
        for idx in match:
            # Check if this atom is part of more than one ring.
            if ring_info.NumAtomRings(idx) > 1:
                fused = True
                break
        if fused:
            # We found a candidate; it has the six-membered ring ketone, 
            # an adjacent double bond and is in a fused ring system.
            return True, ("Molecule contains a six-membered ring ketone (3-oxo group) directly bonded "
                          "to a carbon participating in a C=C double bond (Δ(1)) within a fused ring system.")
    
    return False, "Candidate ketone with adjacent double bond was found but not in a fused six-membered ring system"

# For testing purposes (you can remove or comment out the test cases)
if __name__ == "__main__":
    test_smiles = [
        # True positive (estra-1,5(10)-diene-3,4,17-trione)
        "C[C@]12CC[C@H]3[C@@H](CCC4=C3C=CC(=O)C4=O)[C@@H]1CCC2=O",
        # A false positive example previously misclassified: Isoglabrolide
        "O1C23C(CCC4(C5(C(C6(C(CC5)C(C(O)CC6)(C)C)C)C(=O)C=C42)C)C)(CCC(C3)(C1=O)C)C",
        # A false negative example (paraminabeolide B) – expecting True
        "C[C@@H]([C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3CC[C@]12COC(C)=O"
    ]
    for smi in test_smiles:
        res, reason = is_3_oxo_Delta_1__steroid(smi)
        print("SMILES:", smi)
        print("Classification:", res)
        print("Explanation:", reason)
        print("----------")