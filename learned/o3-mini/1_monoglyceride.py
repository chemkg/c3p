"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride
A 1-monoglyceride is defined as a monoacylglycerol in which the acyl substituent is located at
one of the primary (CH2) positions of the glycerol backbone.
This program first filters molecules with phosphorus (common in phospholipids), then looks for
a unique 1-monoglyceride SMARTS match (which encodes an acyl group, via ester, attached to a 
glycerol backbone that still features free OH groups on the other two primary positions).
Finally, it “grows” the acyl chain and checks that the union of the matched motif and the acyl chain 
comprises the entire molecule. Molecules with extra substituents will be rejected.
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a 1-monoglyceride.
    A 1-monoglyceride must have a glycerol backbone (HO–CH2–CHOH–CH2OH) with exactly one acyl
    (C(=O)–R) group attached as an ester at one of its terminal (primary, CH2) positions.
    
    The algorithm works as follows:
      1. It first creates a molecule from the SMILES string.
      2. It rejects molecules containing phosphorus (P), often present in phospholipids.
      3. It uses a SMARTS pattern that encodes the following motif:
           C(=O)O[CH2]C([OX2H])[CH2][OX2H]
         In this pattern the ester part is:
           • The acyl carbon (C(=O)) [atom index 0],
           • Its carbonyl oxygen (atom index 1),
           • The ester oxygen that links to the glycerol (atom index 2);
         The glycerol backbone is then encoded as:
           • A CH2 that was esterified (atom index 3),
           • A CH with a free OH (atom index 4 with hydroxyl at index 5),
           • A terminal CH2 with a free OH (atom index 6 with hydroxyl at index 7).
      4. If exactly one match is found, we then “grow” the acyl chain from the acyl carbon by
         taking the neighbor of atom0 (the acyl carbon) that is not part of the matched motif.
      5. Finally, we check that the union of the motif (the 8 atoms in the match) plus the full acyl 
         chain comprises the entire molecule. If extra atoms are present, the molecule is not a simple 
         1-monoglyceride.
    
    Args:
      smiles (str): SMILES string representing the molecule.
    
    Returns:
      bool: True if the molecule is a 1-monoglyceride, False otherwise.
      str: A message describing the outcome.
    """
    # Create molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with phosphorus (common in phospholipids)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            return False, "Molecule contains phosphorus; not a 1-monoglyceride"
    
    # Define a SMARTS pattern for the 1-monoglyceride motif.
    # The SMARTS encodes:
    #   [0] Acyl carbon: C(=O)
    #   [1] Carbonyl oxygen (double-bonded)
    #   [2] Ester oxygen that links to glycerol
    #   [3] Esterified glycerol CH2 (primary carbon)
    #   [4] Central glycerol C (with free OH)
    #   [5] Free hydroxyl on central glycerol C
    #   [6] Terminal glycerol CH2 (primary carbon)
    #   [7] Free hydroxyl on terminal glycerol CH2
    # Note: The pattern is written so that the ester group is attached to a CH2 that is clearly
    # part of a glycerol backbone.
    mono_smarts = "C(=O)O[CH2]C([OX2H])[CH2][OX2H]"
    mono_pattern = Chem.MolFromSmarts(mono_smarts)
    if mono_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Find substructure matches (ignoring chirality)
    matches = mol.GetSubstructMatches(mono_pattern, useChirality=False)
    if len(matches) == 0:
        return False, "1-monoglyceride motif not found (ester not attached to a glycerol CH2 with free OH groups)"
    if len(matches) > 1:
        return False, f"Multiple potential 1-monoglyceride motifs found ({len(matches)} matches) in molecule"
    
    # We expect a unique match; get it (a tuple of atom indices)
    match = matches[0]
    # We expect the match to be 8 atoms long (indices 0 through 7 as described in the SMARTS)
    if len(match) != 8:
        return False, f"Match length ({len(match)}) unexpected; expected 8 atoms for the motif"
    
    # Now, extract the acyl chain attached to the acyl carbon.
    # In our motif, match[0] is the acyl carbon (C(=O)).
    acyl_carbon = mol.GetAtomWithIdx(match[0])
    # Identify the neighbor that is not the carbonyl oxygen.
    # (The carbonyl oxygen is in match[1].)
    acyl_root = None
    for nbr in acyl_carbon.GetNeighbors():
        if nbr.GetIdx() not in [match[1]]:  # exclude the carbonyl oxygen (double-bonded O)
            acyl_root = nbr.GetIdx()
            break
    if acyl_root is None:
        return False, "Acyl chain not found attached to the acyl carbon"
    
    # Now, perform a depth-first search (DFS) to identify the entire acyl chain fragment.
    # We start at acyl_root and avoid going back to the motif (specifically, we do not allow going to match[0]).
    acyl_fragment = set()
    stack = [acyl_root]
    while stack:
        current = stack.pop()
        if current in acyl_fragment:
            continue
        acyl_fragment.add(current)
        current_atom = mol.GetAtomWithIdx(current)
        for nbr in current_atom.GetNeighbors():
            # Avoid stepping back to the acyl carbon in the motif.
            if nbr.GetIdx() == match[0]:
                continue
            # Add neighbor if not already part of acyl_fragment.
            if nbr.GetIdx() not in acyl_fragment:
                stack.append(nbr.GetIdx())
    
    # The full motif (the ester + glycerol backbone) is given by the 8 atoms in the match.
    motif_atoms = set(match)
    
    # The union of the motif and the acyl chain is what we expect in a pure 1-monoglyceride.
    union_atoms = motif_atoms.union(acyl_fragment)
    
    all_atoms = set(range(mol.GetNumAtoms()))
    if union_atoms != all_atoms:
        return False, "Extra substituents found outside glycerol monoester motif"
    
    return True, "Found 1-monoglyceride motif: one ester group attached at a primary position of a glycerol backbone"

# Example usage:
if __name__ == "__main__":
    # Test one example: rac-1-monodecanoylglycerol
    test_smiles = "CCCCCCCCCC(=O)OCC(O)CO"
    result, reason = is_1_monoglyceride(test_smiles)
    print(result, reason)