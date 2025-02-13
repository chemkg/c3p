"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: Polar Amino Acids (Free amino acids whose side chain is capable of hydrogen bonding)
Definition: A free amino acid should have a free carboxylic acid group and a free amine group
 (i.e. no peptide bonds) and its side chain (i.e. all atoms connected to the α–carbon except the backbone)
 must contain at least one polar atom (N, O or S). Glycine (side chain = hydrogen) is not considered polar.
"""
from rdkit import Chem
from collections import deque

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a free polar amino acid based on its SMILES string.
    A free polar amino acid is classified as an amino acid that (1) has the expected amino acid
    backbone (free –NH2/–NH3+ and free –COOH (or –COO-)) without any additional peptide (amide)
    bonds and (2) its side chain (i.e. substituents on the α–carbon) contains at least one polar
    atom (a heteroatom: nitrogen, oxygen, or sulfur) capable of hydrogen bonding.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a free polar amino acid, False otherwise
        str: Detailed reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Check for a single free carboxyl group.
    # We use a SMARTS pattern for a carboxylic acid group allowing for either protonated (–OH) or deprotonated ([O-]).
    acid_pattern = Chem.MolFromSmarts("C(=O)[O,OX1-]")  
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Expected exactly one carboxyl group; found {len(acid_matches)}. This may indicate a peptide or non–standard amino acid."

    # (2) Check for peptide bonds (amide bonds). A free amino acid should NOT have a C(=O)–N bond.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Molecule appears to have amide (peptide) bonds and is likely part of a peptide"

    # (3) Identify the amino acid backbone.
    # We allow for stereochemical variations.
    backbone_smarts = ["N[C@@H](*)C(=O)O", "N[C@H](*)C(=O)O"]
    backbone_match = None
    for patt in backbone_smarts:
        patt_mol = Chem.MolFromSmarts(patt)
        matches = mol.GetSubstructMatches(patt_mol)
        if matches:
            # If multiple matches are found, this is ambiguous.
            if len(matches) > 1:
                return False, "Ambiguous backbone matches, molecule may be a peptide"
            backbone_match = matches[0]
            break
    if backbone_match is None:
        return False, "Molecule does not match the expected amino acid backbone pattern"

    # In the matched backbone we assume:
    # backbone_match[0] : amino nitrogen (–NH2 or –NH3+)
    # backbone_match[1] : alpha carbon (C★)
    # backbone_match[2] : carboxyl carbon (C in –COOH)
    # (the pattern "C(=O)O" takes care of both oxygens)
    amino_idx = backbone_match[0]
    alpha_idx = backbone_match[1]
    carboxyl_idx = backbone_match[2]

    # (4) Identify the side chain atoms.
    # Get the alpha carbon atom; its neighbors will include the amino group, the carboxyl carbon, and the side chain.
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    sidechain_start = []
    for neighbor in alpha_atom.GetNeighbors():
        if neighbor.GetIdx() not in (amino_idx, carboxyl_idx):
            sidechain_start.append(neighbor.GetIdx())

    # If there are no side chain atoms then this is glycine (side chain = hydrogen),
    # which is not considered polar.
    if not sidechain_start:
        return False, "Glycine side chain is just a hydrogen and is not polar"

    # Use a breadth-first search (BFS) to collect all atoms in the side chain.
    sidechain_atoms = set()
    queue = deque(sidechain_start)
    while queue:
        current_idx = queue.popleft()
        if current_idx in sidechain_atoms:
            continue
        sidechain_atoms.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nbr in current_atom.GetNeighbors():
            # Do not go back to the alpha carbon (backbone)
            if nbr.GetIdx() == alpha_idx:
                continue
            if nbr.GetIdx() not in sidechain_atoms:
                queue.append(nbr.GetIdx())

    # (5) Check the side chain for polar atoms.
    # Consider atoms nitrogen (7), oxygen (8), and sulfur (16) as polar.
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
    # Several examples for a polar amino acid. Here L-serine is a polar amino acid.
    test_smiles = "N[C@@H](CO)C(O)=O"  # L-serine
    result, reason = is_polar_amino_acid(test_smiles)
    print("Result:", result, "| Reason:", reason)