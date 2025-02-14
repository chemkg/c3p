"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcNumAliphaticCarbocycles

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of a carboxy group of
    a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for ester linkage (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Get the atom indices of all ester linkages found
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Verify that each ester group connects two sufficiently long carbon chains
    for match in ester_matches:
        atom_idx_c1 = match[0]  # Carbon of the carbonyl group
        atom_idx_o = match[2]   # Oxygen in -C(=O)O-
        
        # Check chain connected to atom_idx_c1
        atom_neighbors = [neigh.GetIdx() for neigh in mol.GetAtomWithIdx(atom_idx_c1).GetNeighbors() if neigh.GetIdx() != match[1]]
        c1_chain_length = sum(1 for _ in Chem.rdmolops.bfsTraverse(mol, startAtomIdx=atom_neighbors[0]) if mol.GetAtomWithIdx(_).GetAtomicNum() == 6)
        
        # Check chain connected to atom_idx_o
        atom_neighbors = [neigh.GetIdx() for neigh in mol.GetAtomWithIdx(atom_idx_o).GetNeighbors()]
        o_chain_length = sum(1 for _ in Chem.rdmolops.bfsTraverse(mol, startAtomIdx=atom_neighbors[0]) if mol.GetAtomWithIdx(_).GetAtomicNum() == 6)

        # Ensure both chains are sufficiently long (e.g., >=8 carbons)
        if c1_chain_length < 8 or o_chain_length < 8:
            return False, f"Chains attached to ester group are too short (lengths found: {c1_chain_length} and {o_chain_length})"

    return True, "Contains an ester linkage with sufficiently long carbon chains typical of wax esters"