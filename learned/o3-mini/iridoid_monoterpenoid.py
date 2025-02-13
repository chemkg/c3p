"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
#!/usr/bin/env python
"""
Classifies: Iridoid Monoterpenoid
Definition: Iridoid monoterpenoids are biosynthesized from isoprene and usually have a cyclopentane ring fused
to a six-membered oxygen heterocycle. (Secoiridoids arise when one bond in the cyclopentane ring is cleaved.)
This program detects an "iridoid core" using SMARTS queries.
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines whether the given molecule (SMILES string) is an iridoid monoterpenoid.
    It does so by attempting to match two SMARTS patterns representing:
      1. A classical cyclopenta[c]pyran core (with an oxygen in the fused six-membered ring).
      2. A secoiridoid variant (lacking the oxygen in the six-membered ring).
    In addition, for each match, the number of carbon atoms is counted to be within the 
    expected range for a monoterpenoid core (roughly 7–9 carbons).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as
                     an iridoid monoterpenoid, False otherwise; the second element is a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for the classical cyclopenta[c]pyran core (5-membered ring fused to a 6-membered ring with one oxygen)
    pattern_iridoid = Chem.MolFromSmarts("C1C2CC(C1)OC2")
    # SMARTS for a secoiridoid variant (missing the oxygen in the 6-membered ring)
    pattern_seco = Chem.MolFromSmarts("C1C2CC(C1)C2")
    
    # Helper function: Given a tuple of atom indices (match), count how many are carbons
    def count_carbons(atom_indices):
        return sum(1 for idx in atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # Try classical iridoid pattern
    if mol.HasSubstructMatch(pattern_iridoid):
        matches = mol.GetSubstructMatches(pattern_iridoid)
        for match in matches:
            nC = count_carbons(match)
            if 7 <= nC <= 9:
                return True, ("Found classical iridoid cyclopenta[c]pyran core with appropriate carbon count (%d carbons)" % nC)
        return True, "Potential iridoid core (oxygen-containing) found, but carbon count is outside expected range"
    
    # If classical pattern is not found, try the secoiridoid variant.
    if mol.HasSubstructMatch(pattern_seco):
        matches = mol.GetSubstructMatches(pattern_seco)
        for match in matches:
            nC = count_carbons(match)
            if 7 <= nC <= 9:
                return True, ("Found secoiridoid core variant with appropriate carbon count (%d carbons)" % nC)
        return True, "Potential secoiridoid core found, but carbon count is outside expected range"
    
    return False, "No fused bicyclic iridoid monoterpenoid core was detected"

# (Optional) Test cases using example SMILES from the prompt:
if __name__ == "__main__":
    test_smiles = [
        "BrC1=C(O)C(=C([C@@H]2[C@@](CC[C@@H]2C(C)C)(CC(=O)O)C)C=C1C)C(=O)O",  # Hamigeran L (may be classic or modified)
        "OC1(C2C(CC1O)C(=COC2OC3OC(C(O)C(O)C3O)COC(=O)/C=C/C4=CC=C(O)C=C4)C(O)=O",  # Lippioside I
        "ClCC1(O)C2C(=CC1OC(=O)CC(C)C)C(=COC2OC(=O)CC(C)C)COC(=O)C",  # Valechlorin
        "O=C(OC)[C@@H]1[C@]2(O)[C@]([C@H](OC(=O)C)C([C@@H]2[C@@H](O)C[C@H]1C)(C)C)(CO)C",  # methyl 7alpha-acetoxydeacetylbotryoloate
    ]
    for smi in test_smiles:
        result, reason = is_iridoid_monoterpenoid(smi)
        print("SMILES:", smi)
        print("Result:", result, "| Reason:", reason)
        print("-" * 80)