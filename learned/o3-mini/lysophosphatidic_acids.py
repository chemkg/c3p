"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids
Definition: Any monoacylglycerol phosphate obtained by hydrolytic removal of one of the two acyl groups of any phosphatidic acid or derivatives therein.
This script uses RDKit to check that:
  1. A phosphate group is present (via "[P](=O)(O)(O)").
  2. A non‐ring three‐carbon “glycerol core” is present (via "[CH2;!r]-[CH;!r]-[CH2;!r]").
  3. Exactly one acyl ester group – here defined by an oxygen bound to a carbonyl (and not already bound to a phosphorus) – is attached to the glycerol core.
To further increase specificity, we require that this acyl–oxygen is “near” a phosphate (within 3–5 bonds), as expected from the glycerol–phosphate connection.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    
    The criteria applied are:
      1. The molecule must be valid.
      2. It must contain a phosphate group defined by [P](=O)(O)(O).
      3. It must contain a non‐ring three–carbon fragment corresponding to a glycerol backbone.
      4. It must contain exactly one acyl ester group (an [O]C(=O)[#6] fragment with the oxygen not bound to P)
         and such that the acyl–oxygen lies within 3–5 bonds of a phosphate P (as expected for a glycerol–phosphate connectivity).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a lysophosphatidic acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a phosphate group:
    phosphate_smarts = "[P](=O)(O)(O)"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found"
    
    # Define a SMARTS for a (linear, non‐ring) glycerol core (three connected carbons).
    # This is intended to capture a HO–CH2–CHOH–CH2–OH motif (without ring systems, to avoid sugars).
    glycerol_core_smarts = "[CH2;!r]-[CH;!r]-[CH2;!r]"
    glycerol_core = Chem.MolFromSmarts(glycerol_core_smarts)
    if not mol.HasSubstructMatch(glycerol_core):
        return False, "Glycerol backbone (non-ring 3-carbon chain) not found"
    
    # Find phosphate atoms from our phosphate pattern
    phos_matches = mol.GetSubstructMatches(phosphate_pattern)
    phos_indices = {match[0] for match in phos_matches}  # first atom in the pattern is P
    
    # Define a SMARTS pattern for an acyl ester group.
    # This pattern looks for an oxygen (which is not directly bound to a P) attached to a carbonyl.
    acyl_smarts = "[O;!$([O]-P)]C(=O)[#6]"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No acyl ester group found (expect exactly one fatty acid chain)"
    
    # From the glycerol matches, collect all atom indices present.
    gly_matches = mol.GetSubstructMatches(glycerol_core_smarts)
    gly_atom_idx_set = set()
    for match in gly_matches:
        gly_atom_idx_set.update(match)
    
    # Now select only those acyl ester matches whose ester oxygen is near a phosphate group.
    # We use a distance cutoff: the shortest bond path (number of bonds) between the acyl oxygen
    # and a phosphate P should be between 3 and 5 (typical for a glycerol phosphate connectivity).
    candidate_acyl_count = 0
    for match in acyl_matches:
        # In our acyl SMARTS, match[0] corresponds to the oxygen.
        o_idx = match[0]
        # Optionally, check that this oxygen is connected (by being in close proximity) to the glycerol core.
        # Here we check if it is a member of any glycerol match. (This is a heuristic.)
        if o_idx not in gly_atom_idx_set:
            continue  # acyl ester not attached to the glycerol chain; skip it
        # Now check the shortest path distance from the oxygen to any phosphate P:
        # (i.e. the acyl oxygen should be separated by a few bonds from the phosphate)
        found_near_phos = False
        for p_idx in phos_indices:
            sp = rdmolops.GetShortestPath(mol, o_idx, p_idx)
            # Number of bonds = (number of atoms in the path - 1)
            n_bonds = len(sp) - 1
            if 3 <= n_bonds <= 5:
                found_near_phos = True
                break
        if found_near_phos:
            candidate_acyl_count += 1

    if candidate_acyl_count == 0:
        return False, "No acyl ester group attached to glycerol–phosphate (within expected bond distance) found"
    if candidate_acyl_count > 1:
        return False, f"Multiple acyl ester groups found ({candidate_acyl_count} groups); expected monoacyl structure"
    
    return True, "Molecule has a phosphate group, a suitable non-ring glycerol backbone, and exactly one acyl ester group appropriately linked to the glycerol–phosphate core"

# Example usage (uncomment to test):
# test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCCCCC)(O)(O)=O"  # 1-docosanoyl-glycero-3-phosphate
# result, reason = is_lysophosphatidic_acids(test_smiles)
# print(result, reason)