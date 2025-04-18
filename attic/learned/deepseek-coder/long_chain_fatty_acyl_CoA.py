"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:17408 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA is a CoA molecule with a long-chain fatty acid (C13 to C22) attached via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][SX2][CX3](=[OX1])")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "No CoA moiety found"

    # Look for thioester bond (S-C=O) and ensure it is part of the CoA moiety
    thioester_pattern = Chem.MolFromSmarts("[SX2][CX3](=[OX1])")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond found"

    # Ensure the thioester bond is part of the CoA moiety
    thioester_in_coa = False
    for thioester_match in thioester_matches:
        for coa_match in coa_matches:
            if thioester_match[0] in coa_match and thioester_match[1] in coa_match:
                thioester_in_coa = True
                break
        if thioester_in_coa:
            break
    if not thioester_in_coa:
        return False, "Thioester bond not part of CoA moiety"

    # Find the carbon chain attached to the thioester bond
    fatty_acid_chain = []
    for thioester_match in thioester_matches:
        sulfur_idx = thioester_match[0]
        carbon_idx = thioester_match[1]
        # Traverse the carbon chain starting from the thioester carbon
        current_atom = mol.GetAtomWithIdx(carbon_idx)
        visited = set()
        stack = [(current_atom, 0)]
        while stack:
            atom, depth = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:
                # Exclude carbons that are part of the CoA moiety
                if atom.GetIdx() not in coa_match:
                    fatty_acid_chain.append(atom.GetIdx())
                # Continue traversing through carbon neighbors
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                        stack.append((neighbor, depth + 1))

    # Count unique carbons in the fatty acid chain
    unique_carbons = set(fatty_acid_chain)
    if len(unique_carbons) < 13 or len(unique_carbons) > 24:  # Adjusted range to 13-24
        return False, f"Fatty acid chain length {len(unique_carbons)} is not within the range of 13 to 24 carbons"

    # Check molecular weight - long-chain fatty acyl-CoA typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for long-chain fatty acyl-CoA"

    return True, "Contains CoA moiety with a long-chain fatty acid (C13 to C24) attached via a thioester bond"