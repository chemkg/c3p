"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: Penicillin class
Definition:
  A substituted penam must have:
    • A beta‐lactam ring (4‐membered cyclic amide, SMARTS "N1C(=O)CC1") as an anchor.
    • Two methyl groups at position 2. We proxy this as an S‐bound carbon substituted with two methyl groups 
      (SMARTS "[S]-C([CH3])([CH3])") near the beta‐lactam.
    • A carboxylate substituent at position 3. For our purposes the desired form is the protonated (neutral) one,
      so we search for "C(=O)O" (and not "C(=O)[O-]").
    • A carboxamido group at position 6 (SMARTS "NC(=O)").
      
For each required group we insist that at least one atom in the match lies within 4 bonds of at least one atom
of the beta‐lactam ring. In addition, only molecules that are “neutral” (formal charge of zero) are accepted,
which helps remove counterions and salt forms.
"""

from rdkit import Chem
import numpy as np

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    
    The recognition method:
      1. Parse the SMILES and choose the largest fragment (to remove counterions).
      2. Check that the overall fragment is neutral.
      3. Identify the beta‐lactam ring via the SMARTS "N1C(=O)CC1".
      4. For each of the three groups (S-bound dimethyl, a protonated carboxylate ("C(=O)O")
         and a carboxamido group ("NC(=O)") – verify that at least one match lies within 4 bonds of 
         at least one atom in the beta‐lactam ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as penicillin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If we have more than one fragment (e.g. salts), select the largest
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Ensure that the selected fragment is net neutral.
    if Chem.GetFormalCharge(mol) != 0:
        return False, "Molecule carries a net charge (likely a salt form)"
    
    # Precompute the full distance matrix (bond distances) once
    try:
        dist_matrix = Chem.GetDistanceMatrix(mol)
    except Exception as e:
        return False, f"Error computing distance matrix: {str(e)}"
    
    # Identify beta-lactam ring as the anchor using SMARTS "N1C(=O)CC1"
    beta_lactam_smarts = "N1C(=O)CC1"
    beta_query = Chem.MolFromSmarts(beta_lactam_smarts)
    beta_matches = mol.GetSubstructMatches(beta_query)
    if not beta_matches:
        return False, "Missing beta-lactam ring (SMARTS: N1C(=O)CC1)"
    beta_anchor = set(beta_matches[0])
    
    # Helper function: check if at least one occurrence of a SMARTS pattern is found within max_dist bonds
    # of at least one atom from the beta-lactam anchor.
    def is_close(query_smarts, max_dist=4):
        qmol = Chem.MolFromSmarts(query_smarts)
        if qmol is None:
            return False, f"Invalid SMARTS: {query_smarts}"
        matches = mol.GetSubstructMatches(qmol)
        if not matches:
            return False, f"Group '{query_smarts}' not found"
        for match in matches:
            for atom_idx in match:
                for anchor_idx in beta_anchor:
                    if int(dist_matrix[atom_idx, anchor_idx]) <= max_dist:
                        return True, ""
        return False, f"Group '{query_smarts}' not found within {max_dist} bonds of beta-lactam ring"
    
    # 1. Check for the beta-lactam S-bound dimethyl substituent.
    #    This is our proxy for having two methyl groups at position 2.
    dimethyl_smarts = "[S]-C([CH3])([CH3])"
    close, msg = is_close(dimethyl_smarts)
    if not close:
        return False, "Missing or too distant S-bound dimethyl substituent (expected at position 2): " + msg
    
    # 2. Check for a carboxylate substituent in its neutral (protonated) form.
    #    We require the pattern "C(=O)O". (This avoids matching deprotonated forms.)
    carboxylate_neutral_smarts = "C(=O)O"
    close, msg = is_close(carboxylate_neutral_smarts)
    if not close:
        return False, "Missing protonated carboxylate substituent (C(=O)O) in proximity to beta-lactam ring: " + msg
        
    # 3. Check for a carboxamido group ("NC(=O)")
    carboxamido_smarts = "NC(=O)"
    close, msg = is_close(carboxamido_smarts)
    if not close:
        return False, "Missing carboxamido substituent (NC(=O)) in proximity to beta-lactam ring: " + msg
    
    return True, "Matches penicillin criteria: Contains beta‐lactam, 2,2-dimethyl, protonated carboxylate and carboxamido substituents"

# Example test (uncomment to try):
# test_smiles = "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)CCC[C@H](N)C(O)=O)"
# print(is_penicillin(test_smiles))