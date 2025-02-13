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
  3. Must have a thioester group (S-C(=O)[#6]); from the carbonyl carbon, the contiguous acyl chain (including that carbon) must contain between 6 and 12 carbon atoms.
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
    # Parse SMILES into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Criterion 1: Check for a CoA head group via the adenine substructure.
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine ring (CoA head group not detected)."
    
    # Criterion 2: Check for exactly 4 deprotonated phosphate oxygens 
    # (oxygen atoms with formal charge -1).
    neg_oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if neg_oxygen_count != 4:
        return False, f"Expected 4 negatively charged oxygens from phosphate groups, found {neg_oxygen_count}."
    
    # Criterion 3: Identify the thioester fragment signifying the acyl-CoA unit.
    # SMARTS pattern: S-C(=O)[#6]
    # Note: This will match four atoms: S, the carbonyl carbon, the oxygen (in the C=O), and the acyl chain's first carbon.
    thioester_pattern = Chem.MolFromSmarts("S-C(=O)[#6]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) fragment found."
    
    # Pick the first matching thioester.
    # Unpack four indices: s_idx, carbonyl C index, oxygen index (ignored), and the first acyl chain carbon index.
    try:
        s_idx, coxo_idx, o_idx, acyl_start_idx = thioester_matches[0]
    except Exception as e:
        return False, f"Thioester matching error: {str(e)}"

    # Count the number of contiguous carbon atoms in the acyl chain.
    # We start with the carbonyl carbon and add carbons reachable from the acyl start atom.
    visited = set()

    def dfs(atom_idx):
        """Perform depth-first search to count connected carbon atoms in the acyl chain.
           Do not traverse back to the carbonyl carbon.
        """
        count = 0
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            return 0
        count += 1
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip atoms that are not carbon or have already been counted.
            if nbr_idx not in visited and nbr.GetAtomicNum() == 6:
                count += dfs(nbr_idx)
        return count

    # Begin count from the carbonyl carbon.
    acyl_chain_carbon_count = 1  # Count the carbonyl carbon.
    visited.add(coxo_idx)  # Prevent DFS from including the carbonyl carbon again.
    acyl_chain_carbon_count += dfs(acyl_start_idx)
    
    # Check if the chain falls within the medium-chain length (6-12 carbons inclusive)
    if not (6 <= acyl_chain_carbon_count <= 12):
        return False, f"Acyl chain has {acyl_chain_carbon_count} carbons; expected between 6 and 12 for medium-chain."
    
    return True, f"Detected medium-chain fatty acyl-CoA(4-) with an acyl chain of {acyl_chain_carbon_count} carbons and a valid CoA head group."

# (Optional) To test the function with one of the given examples, uncomment the following lines:
# if __name__ == "__main__":
#     test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/C=C\\CCCCC)=O)=O)O)(C)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
#     result, reason = is_medium_chain_fatty_acyl_CoA_4__(test_smiles)
#     print(result, reason)