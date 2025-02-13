"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: Polar Amino Acids (Free amino acids whose side chain is capable of hydrogen bonding)
Definition: A free polar amino acid should have an unmodified amino acid backbone
(with a free –NH2/–NH3+ and a free carboxyl group) and any substituent on the α–carbon
(i.e., the side chain) must contain at least one polar heteroatom (N, O, or S). Glycine (which has a hydrogen 
as its side chain) is not considered polar.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import deque

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a free polar amino acid based on its SMILES string.
    It performs the following checks:
      (A) The molecule's molecular weight and heavy atom count are within typical amino acid ranges.
      (B) There is exactly one match for the expected amino acid backbone pattern.
          The backbone is defined as a free amine attached to a central (α–) carbon which is bonded to a carboxyl group.
      (C) The side chain – defined as atoms attached to the α–carbon (other than the backbone amine and carboxyl group)
          – is not empty (i.e. the molecule is not glycine) and contains at least one polar atom (N, O, or S).
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule is a free polar amino acid, False otherwise.
      str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (A) Check if molecular weight and heavy atom count are in a typical amino acid range.
    mw = Descriptors.ExactMolWt(mol)
    if not (75 <= mw <= 250):
        return False, f"Molecular weight {mw:.1f} out of range for a free amino acid"
    
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms > 25:
        return False, f"Too many heavy atoms ({heavy_atoms}); likely not a free amino acid"

    # (B) Identify the amino acid backbone.
    # The backbone is defined as: [NH2,NH3+]-C(-*)(C(=O)[O;H1,O-])
    # This pattern allows for the free amine and free carboxylic acid functional groups 
    # without enforcing chirality tags so that molecules with missing/chiral specifications (e.g. tyrosine) match.
    backbone_smarts = "[NH2,NH3+]-C(-*)(C(=O)[O;H1,O-])"
    backbone_patt = Chem.MolFromSmarts(backbone_smarts)
    backbone_matches = mol.GetSubstructMatches(backbone_patt)
    if not backbone_matches:
        return False, "Molecule does not match the expected amino acid backbone pattern"
    if len(backbone_matches) > 1:
        # More than one match may indicate a peptide or ambiguous structure
        return False, "Ambiguous backbone matches; molecule may be part of a peptide"
    
    # Assume the match corresponds to:
    # match[0]: backbone amine atom,
    # match[1]: the α–carbon,
    # match[2]: the carboxyl carbon.
    match = backbone_matches[0]
    backbone_indices = set(match)
    alpha_idx = match[1]  # α–carbon index

    # (C) Identify the side chain.
    # Side chain: all atoms connected to the α–carbon that are not part of the backbone.
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    sidechain_start = []
    for nbr in alpha_atom.GetNeighbors():
        if nbr.GetIdx() not in backbone_indices:
            sidechain_start.append(nbr.GetIdx())
    
    if not sidechain_start:
        return False, "Glycine side chain (only hydrogen) is not polar"
    
    # Use breadth-first search (BFS) to collect all atoms in the side chain.
    sidechain_atoms = set()
    queue = deque(sidechain_start)
    while queue:
        current_idx = queue.popleft()
        if current_idx in sidechain_atoms:
            continue
        sidechain_atoms.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nbr in current_atom.GetNeighbors():
            # Do not traverse back into the α–carbon to avoid looping back into the backbone.
            if nbr.GetIdx() == alpha_idx:
                continue
            if nbr.GetIdx() not in sidechain_atoms:
                queue.append(nbr.GetIdx())
    
    # (D) Check if the side chain contains at least one polar atom.
    # Here, polar atoms are defined as nitrogen (7), oxygen (8), or sulfur (16).
    for idx in sidechain_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() in (7, 8, 16):
            return True, "Side chain contains polar functional group(s) capable of hydrogen bonding"
    
    return False, "Side chain does not contain polar atoms capable of hydrogen bonding"


# Example usage (for testing purposes)
if __name__ == "__main__":
    # Testing with a few examples:
    test_examples = {
        "D-glutamine": "O=C(O)[C@H](N)CCC(=O)N",
        "L-serine": "N[C@@H](CO)C(O)=O",
        "L-glutamic acid": "N[C@@H](CCC(O)=O)C(O)=O",
        "L-lysine": "NCCCC[C@@H](N)C(O)=O",
        "tyrosine": "NC(Cc1ccc(O)cc1)C(O)=O",
        "glycine (should fail)": "NCC(=O)O"
    }
    for name, smi in test_examples.items():
        result, reason = is_polar_amino_acid(smi)
        print(f"SMILES: {smi}\n NAME: {name}\n Result: {result}\n Reason: {reason}\n{'-'*50}")