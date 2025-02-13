"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine
Definition:
  A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
  That is, a glycerol backbone with three –OH groups where the sn-3 hydroxyl is phosphorylated and then
  esterified with L-serine, and the sn-1 and sn-2 hydroxyls are esterified with acyl chains.
  
This program uses RDKit to check that the given SMILES string has a phosphoserine head group 
and at least two acyl ester substituents.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    It verifies that the molecule contains:
      - A phosphoserine head group (a phosphate ester of L-serine attached via a glycerol fragment)
      - Two acyl ester substituents (as fatty acyl chains of at least 4 carbons) attached to positions sn-1 and sn-2.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-sn-phosphatidyl-L-serine, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define two SMARTS patterns for the phosphoserine head group (allowing for different chirality labels).
    ps_pattern1 = Chem.MolFromSmarts("COP(=O)(O)OC[C@H](N)C(=O)O")
    ps_pattern2 = Chem.MolFromSmarts("COP(=O)(O)OC[C@@H](N)C(=O)O")
    
    ps_matches = mol.GetSubstructMatches(ps_pattern1) or mol.GetSubstructMatches(ps_pattern2)
    if not ps_matches:
        return False, "Phosphoserine head group not found"

    # For later exclusion we take the atoms involved in the first phosphoserine match
    ps_atom_indices = set(ps_matches[0])
    
    # Define an ester SMARTS pattern that targets the carbonyl O-ester bond.
    # This pattern [CX3](=O)O matches any ester function.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # To assess acyl chain length, we traverse from the carbonyl carbon.
    def max_chain_length(atom_idx, parent_idx):
        """Recursively calculates the maximum linear chain length (number of C atoms) starting at atom_idx.
           Only traverses bonds to carbon atoms and avoids going back to parent_idx.
        """
        atom = mol.GetAtomWithIdx(atom_idx)
        # Only count carbon atoms
        if atom.GetAtomicNum() != 6:
            return 0
        max_len = 1
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == parent_idx:
                continue
            if nb.GetAtomicNum() == 6:  # if neighbor is carbon, continue along that chain
                current = 1 + max_chain_length(nb_idx, atom_idx)
                if current > max_len:
                    max_len = current
        return max_len

    valid_acyl_count = 0
    # To avoid double-counting the same ester functionality, we record the carbonyl indices seen.
    seen_carbonyls = set()
    
    for match in ester_matches:
        # match[0] is the carbonyl carbon; match[1] is the ester oxygen.
        carbonyl_idx = match[0]
        oxygen_idx = match[1]
        # Exclude if any atom in this ester is part of the phosphoserine head group.
        if (carbonyl_idx in ps_atom_indices) or (oxygen_idx in ps_atom_indices):
            continue

        # In an ester bond R-C(=O)-O-R', the acyl chain is the R group attached to the carbonyl carbon,
        # not the oxygen (which is attached to the glycerol). Identify the neighbor of the carbonyl carbon
        # that is not the oxygen in our match.
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        acyl_start = None
        for nb in carbonyl_atom.GetNeighbors():
            if nb.GetIdx() == oxygen_idx:  # skip the ester oxygen
                continue
            if nb.GetAtomicNum() == 6:  # typically the acyl chain is a carbon
                acyl_start = nb.GetIdx()
                break
        if acyl_start is None:
            continue  # no acyl chain found!
        # Check if this carbonyl was already used in another match.
        if carbonyl_idx in seen_carbonyls:
            continue
        seen_carbonyls.add(carbonyl_idx)

        # Calculate the maximum contiguous chain length starting from the acyl_start atom.
        chain_length = max_chain_length(acyl_start, carbonyl_idx)
        # We require at least 4 carbon atoms in the chain to consider it a fatty acyl substituent.
        if chain_length >= 4:
            valid_acyl_count += 1

    if valid_acyl_count < 2:
        return False, f"Expected at least 2 acyl substituents at the glycerol sn-1 and sn-2 positions; found {valid_acyl_count}"
    
    # As a rough filter, many phosphatidylserines are >500 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low to be a phosphatidylserine"
    
    return True, "Molecule contains a phosphoserine head group with two acyl ester substituents on a glycerol backbone"

# (Optional testing lines)
# test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCC"
# print(is_3_sn_phosphatidyl_L_serine(test_smiles))