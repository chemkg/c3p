"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: Short-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation
of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    
    The algorithm checks:
      - The presence of a thioester group (carbonyl linked to sulfur).
      - That the acyl (fatty acid) fragment attached to the carbonyl carbon is short (no more than 6 carbons).
      - The presence of a CoA moiety by detecting an adenine substructure using several SMARTS patterns
        while ignoring stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES input
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for thioester functionality: carbonyl (C=O) linked to a sulfur atom
    thioester_smarts = "[C](=O)[S]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) bond detected, hence not an acyl-CoA"
    
    # We take the first thioester match.
    # The match tuple corresponds to (acyl_carbon, carbonyl oxygen, sulfur)
    acyl_carbon_idx, _, sulfur_idx = thioester_matches[0]
    
    # --- Count the number of carbon atoms in the acyl (fatty acid) fragment
    # We count carbons connected to the acyl carbon while not crossing into the CoA fragment.
    def dfs_count_atoms(atom_idx, visited):
        count = 0
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            count += 1
        visited.add(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Do not traverse the bond joining to the sulfur (which leads to CoA)
            if nbr_idx == sulfur_idx:
                continue
            if nbr_idx not in visited and nbr.GetAtomicNum() == 6:
                count += dfs_count_atoms(nbr_idx, visited)
        return count

    visited = set()
    acyl_carbons = dfs_count_atoms(acyl_carbon_idx, visited)
    
    # Allowed short-chain fatty acid fragments: acetyl (2C), propionyl (3C), butyryl (4C), valeryl (5C), caproyl (6C)
    if acyl_carbons > 6:
        return False, f"Acyl chain contains {acyl_carbons} carbons, too long for a short-chain fatty acid"
    if acyl_carbons < 2:
        return False, "Acyl chain too short (less than 2 carbons) to be a fatty acid"
    
    # --- Check for a fragment typical of Coenzyme A
    # We look for an adenine moiety by testing several SMARTS patterns.
    adenine_smarts_list = [
        "c1ncnc2ncnc12",     # a generic adenine/purine pattern
        "n1cnc2nc[nH]c2n1",   # alternative ordering
        "c1nc2nc[nH]c2n1"     # another variation
    ]
    adenine_found = False
    for smarts in adenine_smarts_list:
        adenine_pattern = Chem.MolFromSmarts(smarts)
        # Use useChirality=False to overlook stereochemical details in the adenine moiety.
        if mol.HasSubstructMatch(adenine_pattern, useChirality=False):
            adenine_found = True
            break

    if not adenine_found:
        return False, "No adenine fragment typical of CoA detected"
    
    # All checks have passed.
    return True, f"Found thioester bond with an acyl chain of {acyl_carbons} carbons and a CoA moiety"

# Optional: simple testing examples.
if __name__ == "__main__":
    # Test with one of the provided examples: (R)-3-hydroxypentanoyl-CoA.
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C[C@@H](CC)O)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_short_chain_fatty_acyl_CoA(test_smiles)
    print("Result:", result)
    print("Reason:", reason)