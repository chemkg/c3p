"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
Definition: "An acyl-CoA oxoanion that results from deprotonation of the phosphate and diphosphate groups of any medium-chain fatty acyl-CoA; major species at pH 7.3."
Heuristic criteria:
  1. Must have a CoA moiety inferred by the presence of an adenine ring (SMARTS: "c1ncnc2ncnc12").
  2. Must have exactly four negatively charged oxygens (with formal charge -1) that are bonded to a phosphorus atom.
  3. Must have a thioester group (SMARTS: "S-C(=O)[#6]"). Starting from the carbonyl carbon, the longest continuous acyl chain (the “backbone”) – including the carbonyl carbon – must have between 6 and 12 carbon atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a given SMILES string corresponds to a medium-chain fatty acyl-CoA(4-).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): Tuple with True and reason if molecule passes the criteria; else False and reason message.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Criterion 1: Check for CoA head group via an adenine substructure.
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine ring (CoA head group not detected)."
    
    # Criterion 2: Count negatively charged oxygen atoms bonded to phosphorus (phosphate oxygens).
    phosphate_oxygens = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
            # Check if any neighbor is phosphorus (atomic num 15)
            if any(nbr.GetAtomicNum() == 15 for nbr in atom.GetNeighbors()):
                phosphate_oxygens += 1
    if phosphate_oxygens != 4:
        return False, f"Expected 4 negatively charged phosphate oxygens, found {phosphate_oxygens}."
    
    # Criterion 3: Identify the thioester fragment.
    # SMARTS: "S-C(=O)[#6]" matches a sulfur connected to a carbonyl carbon (double bonded to an oxygen) 
    # and then bonded to a carbon.
    thioester_pattern = Chem.MolFromSmarts("S-C(=O)[#6]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) fragment found."
    
    # Use the first thioester match.
    # The match typically returns a tuple: (S_index, carbonyl_C_index, carbonyl_O_index, acyl_start_index)
    try:
        s_idx, carbonyl_idx, oxo_idx, acyl_start_idx = thioester_matches[0]
    except Exception as e:
        return False, f"Thioester matching error: {str(e)}"
    
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Identify the acyl chain branch: the neighbor of the carbonyl that is a carbon and corresponds to acyl_start.
    acyl_branch = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() == acyl_start_idx:
            acyl_branch = nbr
            break
    if acyl_branch is None:
        return False, "Acyl chain branch not found from thioester fragment."
    
    # Define a recursive function to calculate the longest continuous acyl chain from a starting carbon atom.
    # A visited set is used to avoid infinite recursion through cycles.
    def longest_chain_from(atom, visited):
        max_length = 0  # length from this atom in terms of number of carbons (not counting the current one)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:  # only carbon atoms
                continue
            if nbr.GetIdx() in visited:
                continue
            # Continue along this branch. Create a new visited set including the neighbor.
            branch_length = 1 + longest_chain_from(nbr, visited | {nbr.GetIdx()})
            if branch_length > max_length:
                max_length = branch_length
        return max_length

    # Include the carbonyl carbon plus the chain starting at acyl_branch.
    # We start visited with the carbonyl and acyl_start to avoid going back.
    backbone_length = 1 + longest_chain_from(acyl_branch, {carbonyl_atom.GetIdx(), acyl_branch.GetIdx()})
    
    if not (6 <= backbone_length <= 12):
        return False, f"Acyl chain has {backbone_length} carbon atoms; expected between 6 and 12 for medium-chain."
    
    return True, f"Detected medium-chain fatty acyl-CoA(4-) with an acyl chain of {backbone_length} carbons and a valid CoA head group."
    
# Optional test cases (can be uncommented if running as a standalone script)
# if __name__ == "__main__":
#     test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)/C=C/C=C\\CCCCC)=O)=O)O)(C)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
#     result, reason = is_medium_chain_fatty_acyl_CoA_4__(test_smiles)
#     print(result, reason)