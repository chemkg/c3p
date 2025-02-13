"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: Polar Amino Acids (Any amino acid whose side chain is capable of forming one or more hydrogen bonds)
"""
from rdkit import Chem
from collections import deque

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid is defined as an amino acid whose side chain contains at least one heteroatom
    (typically N, O, or S) that can participate in hydrogen bonding.
    
    The algorithm first confirms the presence of an amino acid backbone by matching one of the two common
    stereochemical SMILES patterns: "N[C@@H](*)C(=O)O" or "N[C@H](*)C(=O)O". If found,
    the α-carbon (central carbon) is taken from the match. The side chain is identified as the neighbor(s)
    of the α-carbon that do not belong to the backbone (i.e. not the amino group or the carboxyl group).
    A breadth-first search collects the atoms of the side chain subgraph.
    Finally, if a polar atom is found in the side chain (N, O or S), the molecule is classified as a polar amino acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two backbone patterns to account for stereochemistry variations.
    backbone_patterns = [ "N[C@@H](*)C(=O)O", "N[C@H](*)C(=O)O" ]
    backbone_match = None
    for patt in backbone_patterns:
        pattern = Chem.MolFromSmarts(patt)
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            backbone_match = matches[0]
            break
    if backbone_match is None:
        return False, "Molecule does not match the amino acid backbone pattern"
    
    # In the backbone match, assume the order is:
    # index0: amino nitrogen, index1: alpha carbon, index2: carboxyl carbon, index3: carboxyl oxygen.
    alpha_idx = backbone_match[1]
    amino_idx = backbone_match[0]
    carboxyl_idx = backbone_match[2]
    
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    
    # Identify the side chain neighbor(s): among the α-carbon's neighbors, exclude the backbone atoms.
    sidechain_start = []
    for neighbor in alpha_atom.GetNeighbors():
        if neighbor.GetIdx() not in (amino_idx, carboxyl_idx):
            sidechain_start.append(neighbor.GetIdx())
    
    # If there are no side chain neighbors, the amino acid is glycine (side chain = H) which is not polar.
    if not sidechain_start:
        return False, "Glycine side chain is just a hydrogen and is not polar"
    
    # Collect all atom indices in the side chain using breadth-first search,
    # starting from the side chain neighbor(s), and do not go back to the α-carbon.
    sidechain_atoms = set()
    queue = deque(sidechain_start)
    while queue:
        idx = queue.popleft()
        if idx in sidechain_atoms:
            continue
        sidechain_atoms.add(idx)
        current_atom = mol.GetAtomWithIdx(idx)
        for nbr in current_atom.GetNeighbors():
            # avoid going into the backbone (specifically the alpha carbon)
            if nbr.GetIdx() == alpha_idx:
                continue
            if nbr.GetIdx() not in sidechain_atoms:
                queue.append(nbr.GetIdx())
                
    # Now scan sidechain for polar atoms capable of hydrogen bonding.
    # We include nitrogen (7), oxygen (8) and sulfur (16) as candidates.
    polar_found = False
    for idx in sidechain_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() in (7, 8, 16):
            polar_found = True
            break

    if polar_found:
        return True, "Side chain contains polar functional group(s) capable of hydrogen bonding"
    else:
        return False, "Side chain does not contain polar atoms capable of hydrogen bonding"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test examples: L-serine is a polar amino acid (SMILES: N[C@@H](CO)C(O)=O)
    test_smiles = "N[C@@H](CO)C(O)=O"
    result, reason = is_polar_amino_acid(test_smiles)
    print("Result:", result, "| Reason:", reason)