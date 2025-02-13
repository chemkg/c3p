"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Trienoic fatty acid
Definition: Any polyunsaturated fatty acid (with only C, H, O) that contains exactly three pure C=C bonds
and has a free terminal carboxylic acid group. Additionally, the alkyl chain (starting at the alpha carbon)
must be sufficiently long (at least 8 carbon atoms).
Heuristics used:
  1. Molecule must contain only C, H, and O atoms.
  2. Exactly one free (non‐ester) carboxylic acid group must be present;
     its acid carbon must be terminal (attached to only one carbon).
  3. The alpha carbon (the carbon attached to the acid carbon) must have at least one hydrogen,
     allowing for mild oxidation on that position.
  4. There must be exactly three pure C=C bonds (only between carbons, not counting C=O bonds).
  5. The main acyl chain (starting at the alpha carbon) must be long enough (at least 8 carbons).
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid (by our extended heuristics) must:
      • contain only C, H, and O atoms,
      • possess exactly one free, terminal carboxylic acid group,
      • have an alpha carbon (adjacent to the acid carbon) that is not overly substituted (at least one H),
      • have exactly three pure carbon–carbon double bonds (C=C), and
      • include a sufficiently long carbon chain (at least 8 carbons beyond the acid group).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our criteria, False otherwise.
        str: A reason message explaining the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check that only atoms allowed (C, H, O) are present.
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Contains atoms other than C, H, and O; not a fatty acid"

    # 2. Find carboxylic acid group(s). We look for both protonated and deprotonated forms.
    acid_smarts_list = [
        "[CX3](=O)[OX2H1]",  # e.g. -C(=O)OH
        "[CX3](=O)[O-]"      # e.g. -C(=O)[O-]
    ]
    acid_matches = []
    for smarts in acid_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(patt)
        if matches:
            acid_matches.extend(matches)
    # Require exactly one free acid group.
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxyl groups; expected exactly one free terminal acid group"
    
    acid_match = acid_matches[0]
    # In our SMARTS, the first atom is the carboxyl carbon.
    acidC_idx = acid_match[0]
    acidC = mol.GetAtomWithIdx(acidC_idx)
    
    # 2a. Verify that the acid carbon is terminal: it should have exactly one carbon neighbor.
    carbon_neighbors = [nbr for nbr in acidC.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Acid group is not terminal (acid carbon is attached to more than one carbon)"
    
    # The alpha carbon (neighboring carbon to acid carbon) is expected to be a typical part of the alkyl chain.
    alphaC = carbon_neighbors[0]
    # 2b. Relaxed check: alpha carbon should have at least one hydrogen.
    if alphaC.GetTotalNumHs() < 1:
        return False, "Alpha carbon adjacent to acid group is over/substituted; not a typical fatty acid backbone"
    
    # 3. Count the number of pure carbon–carbon double bonds.
    c_c_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Count only if both atoms are carbon (ignoring C=O bonds).
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                c_c_double_bonds += 1
    if c_c_double_bonds != 3:
        return False, f"Found {c_c_double_bonds} carbon–carbon double bonds; need exactly 3"
    
    # 4. Check that the alkyl chain is sufficiently long.
    # We build an adjacency dictionary for carbon atoms.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_adj = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        a_idx = bond.GetBeginAtomIdx()
        b_idx = bond.GetEndAtomIdx()
        if a_idx in carbon_adj and b_idx in carbon_adj:
            carbon_adj[a_idx].append(b_idx)
            carbon_adj[b_idx].append(a_idx)
    
    # We perform a DFS starting at the alpha carbon
    def dfs(current, visited):
        max_length = 1  # count current carbon
        for neighbor in carbon_adj.get(current, []):
            if neighbor not in visited:
                length = 1 + dfs(neighbor, visited | {neighbor})
                if length > max_length:
                    max_length = length
        return max_length

    alpha_idx = alphaC.GetIdx()
    longest_chain = dfs(alpha_idx, {alpha_idx})
    # Here we require that the chain from the acid group (via the alpha carbon) is long enough.
    if longest_chain < 8:
        return False, f"Longest carbon chain starting from the acid group has only {longest_chain} carbons; chain too short for a fatty acid"

    # All tests passed.
    return True, ("Molecule is a trienoic fatty acid (contains only C, H, and O, exactly one free terminal carboxylic acid group, " +
                  "an alpha carbon that is not over-substituted, exactly three C=C bonds, and a sufficiently long alkyl chain)")

# Example testing; for production these tests can be removed or commented out.
if __name__ == "__main__":
    test_smiles = [
        "OC(=O)CCCCCCC/C=C\\CC(=O)/C=C(\\O)/C=C\\CCO",    # True positive
        "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O",  # True positive
        "CCCC\\C=C/C\\C=C/CC(O)C(O)C\\C=C/CCCC([O-])=O",    # False positive (multiple acid groups or extra substitution)
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCC[C@@H](OO)C(O)=O"      # Example that may be borderline (alpha carbon substitution)
    ]
    for smi in test_smiles:
        result, reason = is_trienoic_fatty_acid(smi)
        print(f"SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")