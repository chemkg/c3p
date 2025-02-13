"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: Naturally occurring prostaglandins (derived from C20 prostanoic acid)

A prostaglandin is defined for this heuristic as a molecule that:
  1. Contains at least 20 carbon atoms.
  2. Contains a non‐aromatic cyclopentane (5‐membered) ring.
  3. Contains a carboxylic acid group (as expected from prostanoic acid or its derivatives).
  4. Contains a long aliphatic chain (at least 4 connected sp3 carbon atoms).
  
Note: This heuristic does not capture every nuance of prostaglandin structure and may have both
false positives and false negatives.
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin (or prostaglandin derivative) based
    on its SMILES string, using the following heuristic criteria:
      1. At least 20 carbon atoms.
      2. A non‐aromatic cyclopentane ring (the prostaglandin cyclopentane core).
      3. Presence of a carboxylic acid group, as expected from prostanoic acid.
      4. Presence of a long aliphatic chain (at least 4 connected sp3 carbon atoms).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets the prostaglandin criteria, False otherwise.
        str: A reason message for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Count carbon atoms and ensure at least 20 are present
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Insufficient carbons: found {carbon_count} but need at least 20 for a prostaglandin"
    
    # Criterion 2: Look for a non-aromatic cyclopentane ring (5-membered ring)
    rings = mol.GetRingInfo().AtomRings()
    has_cyclopentane = False
    for ring in rings:
        if len(ring) == 5:
            # Check that none of the atoms in the ring is aromatic
            if all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                has_cyclopentane = True
                break
    if not has_cyclopentane:
        return False, "No non-aromatic cyclopentane (5-membered) ring found, which is required for a prostaglandin scaffold"
    
    # Criterion 3: Check for the presence of a carboxylic acid group
    carboxylic_smarts = "[CX3](=O)[OX2H]"  # matches -C(=O)OH groups
    carboxylic_pattern = Chem.MolFromSmarts(carboxylic_smarts)
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found, which is required for prostanoic acid derivatives"
    
    # Criterion 4: Look for a long aliphatic side chain.
    # Here we use a SMARTS pattern for at least 4 connected sp3-hybridized carbons.
    chain_smarts = "[CX4]-[CX4]-[CX4]-[CX4]"
    chain_pattern = Chem.MolFromSmarts(chain_smarts)
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficient aliphatic chain (at least 4 connected sp3 carbons) detected"
    
    return True, ("Molecule meets prostaglandin criteria: "
                  "at least 20 carbons, a non-aromatic cyclopentane ring, "
                  "a carboxylic acid group, and a long aliphatic chain.")

# Example usage:
# Uncomment the lines below to test a prostaglandin example
# smiles_example = "CCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1CCCCCCC(O)=O"  # prostaglandin E1
# result, reason = is_prostaglandin(smiles_example)
# print(result, reason)