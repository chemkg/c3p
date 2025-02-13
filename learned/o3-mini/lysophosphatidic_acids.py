"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids
Definition: Any monoacylglycerol phosphate obtained by hydrolytic removal of one of the two acyl groups of any phosphatidic acid or derivatives therein.
This classifier uses RDKit to check that:
  1. A phosphate group is present (via "[P](=O)(O)(O)").
  2. A non‐ring three‐carbon glycerol backbone is present (via "[CH2;!r]-[CH;!r]-[CH2;!r]").
  3. Exactly one acyl ester group (defined via "[O;!$([O]-P)]C(=O)[#6]") is present and appropriately connected (with the acyl oxygen within 3-5 bonds of a phosphate).
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
      4. It must contain exactly one acyl ester group ([O]C(=O)[#6] with the oxygen not bound to P)
         and this acyl ester oxygen lies within 3–5 bonds of a phosphate (supporting the glycerol–phosphate connectivity).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a lysophosphatidic acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a phosphate group.
    phosphate_smarts = "[P](=O)(O)(O)"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found"
    
    # Define the SMARTS pattern for a non‐ring glycerol backbone.
    glycerol_core_smarts = "[CH2;!r]-[CH;!r]-[CH2;!r]"
    glycerol_core = Chem.MolFromSmarts(glycerol_core_smarts)
    if not mol.HasSubstructMatch(glycerol_core):
        return False, "Glycerol backbone (non-ring three-carbon chain) not found"
    
    # Get matches for the glycerol core and collect involved atom indices.
    gly_matches = mol.GetSubstructMatches(glycerol_core)
    gly_atom_idx_set = set()
    for match in gly_matches:
        gly_atom_idx_set.update(match)
    
    # Get the indices of the phosphorus atoms from the phosphate pattern.
    phos_matches = mol.GetSubstructMatches(phosphate_pattern)
    phos_indices = {match[0] for match in phos_matches}  # assume first atom corresponds to P
    
    # Define the SMARTS pattern for the acyl ester group.
    acyl_smarts = "[O;!$([O]-P)]C(=O)[#6]"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    
    if not acyl_matches:
        return False, "No acyl ester group found (expected exactly one acyl chain)"
    
    # Filter the acyl ester matches: the ester oxygen should be part of the glycerol core and near a phosphate.
    candidate_acyl_count = 0
    for match in acyl_matches:
        # In this SMARTS, match[0] is the ester oxygen.
        o_idx = match[0]
        # Check that the oxygen is connected to the glycerol backbone.
        if o_idx not in gly_atom_idx_set:
            continue
        # Check the bond distance from this oxygen to any phosphate phosphorus.
        found_near_phos = False
        for p_idx in phos_indices:
            sp = rdmolops.GetShortestPath(mol, o_idx, p_idx)
            if sp:
                n_bonds = len(sp) - 1  # convert atom count in path to bonds count
                if 3 <= n_bonds <= 5:
                    found_near_phos = True
                    break
        if found_near_phos:
            candidate_acyl_count += 1

    if candidate_acyl_count == 0:
        return False, "No acyl ester group attached appropriately to the glycerol–phosphate motif found"
    if candidate_acyl_count > 1:
        return False, f"Multiple acyl ester groups found ({candidate_acyl_count}); expected only one to match a lysophosphatidic acid structure"
    
    return True, "Molecule has a phosphate group, a suitable non-ring glycerol backbone, and exactly one acyl ester group correctly attached to form a lysophosphatidic acid"

# Example usage (uncomment the lines below to test):
# test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCCCCC)(O)(O)=O"  # Example: 1-docosanoyl-glycero-3-phosphate
# result, explanation = is_lysophosphatidic_acids(test_smiles)
# print(result, explanation)