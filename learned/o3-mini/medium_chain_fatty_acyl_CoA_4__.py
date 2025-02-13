"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: medium-chain fatty acyl-CoA(4-)
Definition: "An acyl-CoA oxoanion that results from deprotonation of the phosphate and diphosphate groups of any medium-chain fatty acyl-CoA; major species at pH 7.3."
Heuristic criteria:
  1. Must have a CoA moiety inferred by the presence of an adenine ring (SMARTS: "c1ncnc2ncnc12").
  2. Must have exactly four negatively charged oxygen atoms (formal charge –1) from the phosphate/diphosphate groups.
  3. Must have a thioester (S–C(=O)[#6]) group; from the carbonyl carbon, the contiguous acyl chain (counting it as well as the “alpha carbon”) must contain between 6 and 12 carbon atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if the given SMILES string corresponds to a medium-chain fatty acyl-CoA(4-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a medium-chain fatty acyl-CoA(4-), False otherwise.
        str: A reason describing the outcome.
    """
    # Parse SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Criterion 1: Check for a CoA head group via the adenine substructure.
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine ring (CoA head group not detected)."

    # Criterion 2: Check for 4 deprotonated phosphate oxygens (atoms with formal charge -1 on O).
    neg_oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if neg_oxygen_count != 4:
        return False, f"Expected 4 negatively charged oxygens from phosphate groups, found {neg_oxygen_count}."

    # Criterion 3: Look for the thioester group that defines the acyl chain.
    # SMARTS pattern: sulfur attached to a carbonyl carbon which is attached to a carbon (the acyl chain)
    thioester_smarts = "S-C(=O)[#6]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) fragment found."

    # We will choose the first matching thioester.
    # In the SMARTS, match indices: 0 -> S, 1 -> carbonyl C, 2 -> first acyl chain carbon.
    s_idx, coxo_idx, acyl_start_idx = thioester_matches[0]

    # To count the carbons in the fatty acyl chain we include the carbonyl carbon plus all contiguous carbons
    # that extend from the acyl_start atom but do not go back into the CoA moiety.
    visited = set()
    def dfs(atom_idx):
        """Depth-first search to count contiguous carbon atoms (excluding bonds back to the carbonyl)."""
        count = 0
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        # Only follow carbon atoms (atomic number 6)
        if atom.GetAtomicNum() != 6:
            return 0
        count += 1
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Ensure we only follow carbons and do not go back to the carbonyl carbon (which is already counted)
            if nbr_idx not in visited and nbr.GetAtomicNum() == 6:
                count += dfs(nbr_idx)
        return count

    # Start the acyl chain count from the acyl_start atom.
    # Also include the carbonyl carbon (which is part of the fatty acyl chain).
    acyl_chain_carbon_count = 1  # for the carbonyl carbon
    visited.add(coxo_idx)  # do not go back to the carbonyl carbon during DFS
    acyl_chain_carbon_count += dfs(acyl_start_idx)

    # For medium-chain acyl groups, we require between 6 and 12 carbons (inclusive)
    if not (6 <= acyl_chain_carbon_count <= 12):
        return False, f"Acyl chain contains {acyl_chain_carbon_count} carbons; expected between 6 and 12 for medium-chain."

    return True, f"Detected medium-chain fatty acyl-CoA(4-) with an acyl chain of {acyl_chain_carbon_count} carbons and CoA head group."

# (Optional) For testing purposes you could uncomment and test one of the examples:
# if __name__ == "__main__":
#     test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/C=C\\CCCCC)=O)=O)O)(C)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
#     result, reason = is_medium_chain_fatty_acyl_CoA_4__(test_smiles)
#     print(result, reason)