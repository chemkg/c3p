"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine
Definition: An sn-glycerophosphoserine compound having an acyl substituent at the 1-hydroxy position.
The algorithm requires:
  1. The presence of a serine headgroup fragment (approximated as an alcohol attached to a CH(N)C(=O)O unit).
  2. The presence of an acylated glycerol fragment (approximated as an alcohol linked to a -CH2CH(OH)COC(=O)[#6] group).
  3. That the oxygen that bridges these fragments is attached to the same phosphorus atom.
  4. That the acyl chain (extending from the carbonyl carbon) is “long” enough (we require at least 4 contiguous carbons).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def count_contiguous_carbon_chain(mol, start_idx, exclude_idx):
    """
    Simple DFS to count the number of carbon atoms connected to start_idx,
    excluding a given neighbor (for instance the carbonyl carbon) to avoid
    backtracking. Returns the maximum chain length.
    """
    visited = set()
    def dfs(atom_idx, parent):
        visited.add(atom_idx)
        max_length = 1  # count current carbon
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nid = nbr.GetIdx()
            if nid == parent or nid == exclude_idx or nid in visited:
                continue
            if nbr.GetAtomicNum() == 6:  # carbon
                length = 1 + dfs(nid, atom_idx)
                if length > max_length:
                    max_length = length
        return max_length
    return dfs(start_idx, exclude_idx)

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    
    The procedure is:
      - Look for exactly one serine headgroup fragment.
         (SMARTS pattern: OCC(N)C(=O)O)
      - Look for exactly one acylated glycerol fragment at sn-1.
         (SMARTS pattern: OCC(O)COC(=O)[#6])
      - For each, obtain the bridging oxygen (assumed to be the first atom in the pattern) and verify that
        it is directly attached to a phosphorus atom.
      - Ensure that the same phosphorus atom is bridging both fragments.
      - Additionally, follow the acyl chain (from the carbon attached to the carbonyl group)
        and ensure that it contains at least 4 contiguous carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fits the 1-acyl-sn-glycero-3-phosphoserine class, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the two key fragments.
    # We do not include stereochemistry in the SMARTS to allow for variant representations.
    serine_pattern = Chem.MolFromSmarts("OCC(N)C(=O)O")
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)[#6]")
    
    serine_matches = mol.GetSubstructMatches(serine_pattern)
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    
    if not serine_matches:
        return False, "Serine headgroup fragment (OCC(N)C(=O)O) not found"
    if not glycerol_matches:
        return False, "Acylated glycerol fragment (OCC(O)COC(=O)[#6]) not found"
    
    # For clarity we require exactly one match for each fragment.
    if len(serine_matches) > 1:
        return False, f"Ambiguous serine fragment: {len(serine_matches)} matches found"
    if len(glycerol_matches) > 1:
        return False, f"Ambiguous acylated glycerol fragment: {len(glycerol_matches)} matches found"
    
    # For each match, assume the first atom in the match is the bridging oxygen that should be attached to phosphorus.
    serine_bridge_idx = serine_matches[0][0]
    glycerol_bridge_idx = glycerol_matches[0][0]
    
    # Obtain phosphorus neighbors attached to these bridging oxygens.
    serine_p_set = {nbr.GetIdx() for nbr in mol.GetAtomWithIdx(serine_bridge_idx).GetNeighbors() if nbr.GetAtomicNum() == 15}
    if not serine_p_set:
        return False, "Serine fragment found but its bridging oxygen is not attached to any phosphorus atom"
    
    glycerol_p_set = {nbr.GetIdx() for nbr in mol.GetAtomWithIdx(glycerol_bridge_idx).GetNeighbors() if nbr.GetAtomicNum() == 15}
    if not glycerol_p_set:
        return False, "Acylated glycerol fragment found but its bridging oxygen is not attached to any phosphorus atom"
    
    # Check that there is at least one phosphorus atom common to both sets.
    common_p = serine_p_set.intersection(glycerol_p_set)
    if not common_p:
        return False, "Bridging oxygens from the serine and glycerol fragments are not attached to the same phosphorus atom"
    
    # Optionally, check that the acyl chain is of sufficient length.
    # In the glycerol pattern "OCC(O)COC(=O)[#6]", the last atom ([#6]) is the start of the acyl chain.
    # We will attempt to count contiguous carbon atoms beyond that atom.
    glycerol_match = glycerol_matches[0]
    if len(glycerol_match) < 6:
        return False, "Acylated glycerol fragment match is incomplete"
    acyl_start_idx = glycerol_match[5]  # index of the first carbon of the acyl chain
    # We exclude the carbonyl carbon (which is the previous atom in the pattern, index 4).
    carbonyl_idx = glycerol_match[4]
    chain_length = count_contiguous_carbon_chain(mol, acyl_start_idx, carbonyl_idx)
    if chain_length < 4:
        return False, f"Acyl chain too short (chain length {chain_length})"
    
    return True, "Contains a phosphoserine headgroup with an acylated glycerol fragment at the sn-1 position"

# For manual testing:
if __name__ == "__main__":
    # Example correct SMILES:
    test_smiles = [
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",  # PS(17:1(9Z)/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",         # PS(16:0/0:0)
        "P(OCC(O)COC(=O)CCCCCCCCCCCCCCCCC)(OCC(N)C(O)=O)(O)=O",                 # 1-Stearoylglycerophosphoserine
    ]
    for smi in test_smiles:
        result, reason = is_1_acyl_sn_glycero_3_phosphoserine(smi)
        print("SMILES:", smi)
        print("Result:", result, "|", reason)
        print("------------------------------------------------")