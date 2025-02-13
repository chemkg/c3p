"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: Naturally occurring prostaglandins (derived from C20 prostanoic acid)

A prostaglandin is defined for this heuristic as a molecule that:
  1. Contains at least 20 carbon atoms.
  2. Contains a non‐aromatic cyclopentane (5-membered) ring.
  3. Contains a carboxyl functionality (free acid or ester) that is directly attached to that cyclopentane ring.
  4. Contains a long aliphatic chain (at least 4 connected sp3 carbon atoms) as a side‐chain.

Note: This heuristic does not capture every nuance of prostaglandin chemistry;
it has been modified to reduce false positives by requiring that the carboxyl/carboxylate group
be connected to the cyclopentane core, and to accept esterified derivatives.
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin (or prostaglandin derivative) based on its
    SMILES string. It uses the following improved heuristic criteria:
      1. The molecule must contain at least 20 carbon atoms.
      2. It must contain at least one non‐aromatic cyclopentane ring.
      3. It must contain a carboxyl or carboxylate ester group that is directly attached to the cyclopentane ring.
      4. It must contain at least one long aliphatic chain of 4 connected sp3 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the prostaglandin criteria, False otherwise.
        str: A reason message for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Criterion 1: At least 20 carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Insufficient carbons: found {carbon_count} but need at least 20 for a prostaglandin"
    
    # Criterion 2: Identify a non-aromatic cyclopentane (5-membered) ring.
    rings = mol.GetRingInfo().AtomRings()
    cyclopentane_rings = []
    for ring in rings:
        if len(ring) == 5 and all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            cyclopentane_rings.append(ring)
    if not cyclopentane_rings:
        return False, "No non-aromatic cyclopentane (5-membered) ring found, which is required for a prostaglandin scaffold"

    # Criterion 3:
    # Look for a carboxyl / carboxylate ester group.
    # The SMARTS below matches a carbonyl (C(=O)) connected to an oxygen in either a -OH or -OR configuration.
    carboxyl_smarts = "[CX3](=O)O"  # will match -COOH and -COOR fragments.
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl or ester group (-C(=O)O-) found, which is required for prostanoic acid derivatives"
    
    # Now ensure that at least one of the carboxyl groups is directly attached to a cyclopentane ring.
    attached = False
    for ring in cyclopentane_rings:
        for match in carboxyl_matches:
            # Here match[0] is the carbonyl carbon.
            carboxyl_c = match[0]
            atom = mol.GetAtomWithIdx(carboxyl_c)
            # Check all neighbors of the carboxyl carbon:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    attached = True
                    break
            if attached:
                break
        if attached:
            break
    if not attached:
        return False, "Found a carboxyl/ester group, but none is directly attached to the cyclopentane ring"
    
    # Criterion 4: Look for a long aliphatic chain.
    # This SMARTS matches four connected non-aromatic, sp3 carbons.
    chain_smarts = "[CX4]-[CX4]-[CX4]-[CX4]"
    chain_pattern = Chem.MolFromSmarts(chain_smarts)
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficient aliphatic chain (at least 4 connected sp3 carbons) detected"
    
    return True, ("Molecule meets prostaglandin criteria: "
                  "at least 20 carbons, a non-aromatic cyclopentane ring with an attached carboxyl/ester group, "
                  "and a long aliphatic chain.")

# Example usage:
# Uncomment the lines below to test an example prostaglandin SMILES
# smiles_example = "CCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1CCCCCCC(O)=O"  # prostaglandin E1
# result, reason = is_prostaglandin(smiles_example)
# print(result, reason)