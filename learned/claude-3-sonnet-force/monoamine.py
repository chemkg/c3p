"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: CHEBI:35457 monoamine

A monoamine is an aralylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain.
Monoamines are derived from aromatic amino acids like phenylalanine, tyrosine, tryptophan, and the thyroid hormones by the action of aromatic amino acid decarboxylase enzymes.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a monoamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    aromatic_ring_matches = mol.GetSubstructMatches(aromatic_ring_pattern)
    if not aromatic_ring_matches:
        return False, "No aromatic ring found"

    # Look for amino group (-NH2, -NH3+, -NHR, -NR2, -NR3+, =N-, =N+)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino group found"

    # Check if the amino group is connected to an aromatic ring via a two-carbon chain
    for amino_idx in amino_matches[0]:
        for ring_idx in aromatic_ring_matches[0]:
            path = Chem.FindAllPathsOfLengthN(mol, amino_idx, ring_idx, 3, useBonds=True)
            if path:
                chain_atoms = [mol.GetAtomWithIdx(idx) for idx in path[0][1:-1]]
                if len(chain_atoms) == 2 and all(atom.GetAtomicNum() == 6 for atom in chain_atoms):
                    return True, "Monoamine structure detected"

    return False, "Amino group not connected to aromatic ring via two-carbon chain"