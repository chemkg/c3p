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
        (bool, str): True and a reason message if the molecule is classified as a medium-chain fatty acyl-CoA(4-);
                     otherwise False and a corresponding failure message.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Criterion 1: Check for CoA head group via an adenine substructure.
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine ring (CoA head group not detected)."
    
    # Criterion 2: Count negative oxygens that are part of phosphate groups.
    # They are defined as oxygen atoms (atomic num 8) with formal charge -1 that are bonded to a phosphorus atom.
    phosphate_oxygens = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
            # Check if any neighbor is phosphorus (atomic number 15)
            if any(nbr.GetAtomicNum() == 15 for nbr in atom.GetNeighbors()):
                phosphate_oxygens += 1
    if phosphate_oxygens != 4:
        return False, f"Expected 4 negatively charged phosphate oxygens, found {phosphate_oxygens}."
    
    # Criterion 3: Look for the thioester fragment.
    # SMARTS: "S-C(=O)[#6]" matches a sulfur connected to a carbonyl carbon (double bonded to O) which is then bound to an alkyl carbon.
    thioester_pattern = Chem.MolFromSmarts("S-C(=O)[#6]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) fragment found."
    
    # Use the first thioester match.
    # The match returns a tuple of indices corresponding to (S, carbonyl C, carbonyl O, first acyl chain C).
    try:
        s_idx, carbonyl_idx, oxo_idx, acyl_start_idx = thioester_matches[0]
    except Exception as e:
        return False, f"Thioester matching error: {str(e)}"
    
    # We now wish to measure the continuous length of the acyl chain backbone.
    # By definition, the acyl chain backbone starts at the carbonyl carbon (which is part of the thioester group)
    # and continues from the acyl_start atom (the alkyl group attached to the carbonyl C) in the longest path.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    # In the thioester group, carbonyl is bonded to the S atom and to the acyl_start carbon.
    # We want only the acyl chain branch (skip the sulfur branch).
    acyl_branch = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() == acyl_start_idx:
            acyl_branch = nbr
            break
    if acyl_branch is None:
        return False, "Acyl chain branch not found from thioester fragment."
    
    # We now compute the longest continuous (simple) path in the acyl chain starting from acyl_start.
    # We define a recursive function that (given a current carbon atom and the atom index it came from) returns
    # the length (number of carbon atoms) of the longest path from this atom.
    def longest_chain_from(atom, parent_idx):
        max_length = 1  # count this atom
        for nbr in atom.GetNeighbors():
            # Only proceed if neighbor is carbon and is not the parent.
            if nbr.GetIdx() == parent_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            length = 1 + longest_chain_from(nbr, atom.GetIdx())
            if length > max_length:
                max_length = length
        return max_length

    # The overall acyl chain length includes the carbonyl carbon plus the longest chain starting at acyl_start.
    backbone_length = 1 + longest_chain_from(acyl_branch, carbonyl_idx)
    
    # Check if the acyl chain backbone length is in the desired range.
    if not (6 <= backbone_length <= 12):
        return False, f"Acyl chain has {backbone_length} carbons; expected between 6 and 12 for medium-chain."
    
    return True, f"Detected medium-chain fatty acyl-CoA(4-) with an acyl chain of {backbone_length} carbons and a valid CoA head group."

# (Optional testing)
# if __name__ == "__main__":
#     test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/C=C\\CCCCC)=O)=O)O)(C)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
#     result, reason = is_medium_chain_fatty_acyl_CoA_4__(test_smiles)
#     print(result, reason)