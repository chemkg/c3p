"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: Polar Amino Acids (Free amino acids whose side chain is capable of hydrogen bonding)
Definition: A free polar amino acid should have an unmodified amino acid backbone
(with a free –NH2/–NH3+ and a free carboxyl group) and any substituent on the α–carbon
(i.e. the side chain) must contain at least one polar heteroatom (N, O, or S).
Glycine (which has a hydrogen for a side chain) is not considered polar.
Additionally, the molecule’s overall size should be within the typical range for amino acids.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import deque

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a free polar amino acid based on its SMILES string.
    It performs the following checks:
      (1) The molecule must have a molecular weight and heavy atom count consistent with a free amino acid.
      (2) It must have exactly one amino acid backbone pattern match.
          The backbone is defined by N - α–C - C(=O)[O,OX1-],
          allowing for either chirality at the alpha carbon.
      (3) Once the alpha carbon is identified, its neighbors (other than the backbone atoms)
          are taken as the side chain. That side chain must include at least one polar atom (N, O, or S).
          (A side chain with no extra atoms, as in glycine, is not polar.)
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      bool: True if the molecule is a free polar amino acid, False otherwise.
      str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (A) Check if molecular weight and heavy atom count are in a typical amino acid range.
    mw = Descriptors.ExactMolWt(mol)
    if not (75 <= mw <= 250):
        return False, f"Molecular weight {mw:.1f} out of range for free amino acids"
    
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms > 25:
        return False, f"Too many heavy atoms ({heavy_atoms}); likely not a free amino acid"

    # (B) Identify the amino acid backbone.
    # We allow for either chirality at the α–carbon so we use two patterns.
    # The expected backbone: a nitrogen attached to a chiral alpha carbon then a carboxyl group.
    backbone_smarts_list = [
        "N[C@@H](*)C(=O)[O,OX1-]",  # alpha carbon chiral specification (S or R)
        "N[C@H](*)C(=O)[O,OX1-]"
    ]
    
    backbone_match = None
    match_count = 0
    for smarts in backbone_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(patt)
        if matches:
            match_count += len(matches)
            # if more than one match is found, that's ambiguous
            if len(matches) > 1:
                return False, "Ambiguous backbone matches; molecule may be a peptide"
            backbone_match = matches[0]
            
    if backbone_match is None or match_count == 0:
        return False, "Molecule does not match the expected amino acid backbone pattern"
    
    # Backbone match ordering assumption:
    # backbone_match[0] : Amino nitrogen
    # backbone_match[1] : Alpha carbon (C*)
    # backbone_match[2] : Carboxyl carbon (of -COOH)
    # (Other atoms in the COOH group are not used to determine the side chain.)
    amino_idx = backbone_match[0]
    alpha_idx = backbone_match[1]
    backbone_indices = set(backbone_match)  # Include backbone atoms from the match
    
    # (C) Check for peptide (amide) bonds beyond the expected backbone.
    # A free amino acid should not have additional C(=O)N bonds.
    # We look for any amide pattern anywhere in the molecule.
    amide_patt = Chem.MolFromSmarts("C(=O)N")
    # However, if the only amide match corresponds to the backbone pattern, that's acceptable.
    all_amide_matches = mol.GetSubstructMatches(amide_patt)
    if all_amide_matches:
        # Count amide bonds not part of the backbone.
        extra_amide = False
        for match in all_amide_matches:
            # If the carbon in the match is not the backbone carboxyl carbon, then flag as extra.
            if match[0] not in backbone_indices:
                extra_amide = True
                break
        if extra_amide:
            return False, "Molecule appears to have extra amide (peptide) bonds; likely part of a peptide"
    
    # (D) Identify the side chain.
    # The side chain is defined as the atoms (and all connected atoms) that are attached to the α–carbon,
    # except for the backbone atoms (the amino nitrogen and the carboxyl group already identified).
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    sidechain_start = []
    for nbr in alpha_atom.GetNeighbors():
        if nbr.GetIdx() not in backbone_indices:
            sidechain_start.append(nbr.GetIdx())
    
    if not sidechain_start:
        return False, "Glycine side chain (only hydrogen) is not polar"
    
    # Use breadth-first search (BFS) to collect all atoms of the side chain.
    sidechain_atoms = set()
    queue = deque(sidechain_start)
    while queue:
        current = queue.popleft()
        if current in sidechain_atoms:
            continue
        sidechain_atoms.add(current)
        atom = mol.GetAtomWithIdx(current)
        for nbr in atom.GetNeighbors():
            # Do not traverse back into the backbone (alpha atom)
            if nbr.GetIdx() == alpha_idx:
                continue
            if nbr.GetIdx() not in sidechain_atoms:
                queue.append(nbr.GetIdx())

    # (E) Check if the side chain contains any polar atoms.
    # Consider nitrogen (7), oxygen (8), and sulfur (16) as indicators of polarity.
    for idx in sidechain_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() in (7, 8, 16):
            return True, "Side chain contains polar functional group(s) capable of hydrogen bonding"
    
    return False, "Side chain does not contain polar atoms capable of hydrogen bonding"


# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with a known polar amino acid: L-serine
    test_smiles = "N[C@@H](CO)C(O)=O"  # L-serine
    result, reason = is_polar_amino_acid(test_smiles)
    print("Result:", result, "| Reason:", reason)