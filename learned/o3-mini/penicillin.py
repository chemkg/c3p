"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: Penicillin class
Definition:
  Any member of the group of substituted penams containing two methyl substituents at position 2,
  a carboxylate substituent at position 3 and a carboxamido group at position 6.
Heuristic:
  • Use the beta‐lactam ring as an anchor (SMARTS: "N1C(=O)CC1").
  • Check for an S‐bound carbon with two methyl groups (SMARTS: "[S]-C([CH3])([CH3])").
  • Check for a carboxylate substituent (either "C(=O)[O-]" or "C(=O)O").
  • Check for a carboxamido group (SMARTS: "NC(=O)").
  • Each required functionality must be within 4 bonds of at least one atom in the beta‐lactam ring.
The method first selects the largest fragment (to remove counterions) and then computes the full distance matrix
to avoid repeated shortest path calls that can trigger internal errors.
"""

from rdkit import Chem
import numpy as np

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    
    The function identifies the beta-lactam ring (SMARTS "N1C(=O)CC1"),
    and then ensures that within four bonds of at least one of its atoms there exists:
    1. A dimethyl group on an S-bound carbon ("[S]-C([CH3])([CH3])")
    2. A carboxylate substituent (either "C(=O)[O-]" or "C(=O)O")
    3. A carboxamido group ("NC(=O)")
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as penicillin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If the molecule has multiple fragments (e.g. salts), choose the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Compute the complete distance matrix (bond distances) once.
    try:
        dist_matrix = Chem.GetDistanceMatrix(mol)
    except Exception as e:
        return False, f"Error computing distance matrix: {str(e)}"
    
    # Define the beta-lactam ring SMARTS and find matches.
    beta_lactam_smarts = "N1C(=O)CC1"
    beta_query = Chem.MolFromSmarts(beta_lactam_smarts)
    beta_matches = mol.GetSubstructMatches(beta_query)
    if not beta_matches:
        return False, "Missing beta-lactam ring (4-membered cyclic amide)"
    # Use the first match as our anchor.
    beta_anchor = set(beta_matches[0])
    
    # Helper function: For a given SMARTS pattern, check that at least one atom
    # in one of the matching groups is within max_dist bonds (using the precomputed matrix)
    # of at least one of the beta_lactam anchor atoms.
    def is_close(query_smarts, max_dist=4):
        query_mol = Chem.MolFromSmarts(query_smarts)
        if query_mol is None:
            return False, f"Invalid SMARTS: {query_smarts}"
        matches = mol.GetSubstructMatches(query_mol)
        if not matches:
            return False, f"Missing required group: {query_smarts}"
        for match in matches:
            for atom_idx in match:
                for anchor_idx in beta_anchor:
                    if dist_matrix[atom_idx, anchor_idx] <= max_dist:
                        return True, ""
        # If no atom in any match is within the specified bond distance.
        return False, f"Group '{query_smarts}' not found within {max_dist} bonds of the beta-lactam ring"
    
    # 1. S-bound dimethyl substituent (proxy for 2,2-dimethyl) using SMARTS:
    dimethyl_smarts = "[S]-C([CH3])([CH3])"
    close, msg = is_close(dimethyl_smarts)
    if not close:
        return False, "Missing or too distant S-bound dimethyl substituent (expected at position 2): " + msg
    
    # 2. Carboxylate substituent.
    # Allow for protonated ("C(=O)O") or deprotonated ("C(=O)[O-]") forms.
    carboxylate_smarts1 = "C(=O)[O-]"
    carboxylate_smarts2 = "C(=O)O"
    close1, msg1 = is_close(carboxylate_smarts1)
    close2, msg2 = is_close(carboxylate_smarts2)
    if not (close1 or close2):
        return False, "Missing carboxylate substituent (C(=O)O or C(=O)[O-]) in proximity to beta-lactam ring"
        
    # 3. Carboxamido group: "NC(=O)"
    carboxamido_smarts = "NC(=O)"
    close, msg = is_close(carboxamido_smarts)
    if not close:
        return False, "Missing carboxamido group (NC(=O)) in proximity to beta-lactam ring: " + msg
    
    return True, "Matches penicillin criteria: Contains beta‐lactam, 2,2-dimethyl, carboxylate and carboxamido substituents"

# You can test the function with one of the examples. For instance:
# test_smiles = "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)CCC[C@H](N)C(O)=O)"
# print(is_penicillin(test_smiles))