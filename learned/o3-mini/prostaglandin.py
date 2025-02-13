"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: Naturally occurring prostaglandins (derived from C20 prostanoic acid)
A prostaglandin is typically built upon a cyclopentane (5‐membered) ring,
has overall at least 20 carbon atoms, and features a long aliphatic side chain.
Note: This heuristic does not capture all possible structures but is a first‐order filter.
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin (or prostaglandin derivative) based on its SMILES string.
    
    The heuristic criteria are:
      1. The molecule must have at least 20 carbon atoms (reflecting the C20 origin).
      2. The molecule must contain a cyclopentane ring (i.e. at least one ring of exactly 5 atoms).
      3. The molecule must have a long aliphatic chain (we search for a chain of 4 connected aliphatic carbons)
         as one indicator for a prostanoid side chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the prostaglandin criteria, False otherwise.
        str: A reason message for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Insufficient carbons: found {carbon_count} but need at least 20 for a prostaglandin"

    # Check for a cyclopentane ring (a 5-membered ring)
    ring_info = mol.GetRingInfo().AtomRings()
    has_cyclopentane = any(len(ring) == 5 for ring in ring_info)
    if not has_cyclopentane:
        return False, "No cyclopentane (5-membered) ring found, which is required for a prostaglandin scaffold"
    
    # Check for the presence of a long aliphatic chain.
    # We search for a chain of 4 connected sp3 carbon atoms as a heuristic for a side chain.
    chain_pattern = Chem.MolFromSmarts("[CX4]-[CX4]-[CX4]-[CX4]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficient aliphatic chain (at least 4 connected sp3 carbons) detected"
    
    # If all criteria are met, classify as prostaglandin
    return True, ("Molecule has at least 20 carbons, a cyclopentane ring, "
                  "and a long aliphatic chain, consistent with prostaglandin structures.")
                    
# Example usage: Uncomment the following lines to test on one of the provided SMILES.
# smiles_example = "CCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1CCCCCCC(O)=O"  # prostaglandin E1
# result, reason = is_prostaglandin(smiles_example)
# print(result, reason)