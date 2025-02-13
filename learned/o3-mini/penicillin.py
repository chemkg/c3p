"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: Penicillin class
Definition:
  Any member of the group of substituted penams containing two methyl substituents at position 2,
  a carboxylate substituent at position 3 and a carboxamido group at position 6.
Heuristic:
  • Find a beta‐lactam ring, defined here by the SMARTS "N1C(=O)CC1".
  • Locate a sulfur–bound carbon bearing two methyl groups, as a proxy for the C2 dimethyl substituents, via "[S]-C([CH3])([CH3])".
  • Find a carboxylate substituent (accepting both protonated and deprotonated forms).
  • Find a carboxamido group (NC(=O)).
  • In addition, require that each of these substituents is “close” to the beta–lactam ring (within 4 bonds).
Note: The method first removes small fragments (e.g. counterions) and works on the largest fragment.
"""

from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    
    The function uses the beta-lactam ring as an anchor and checks that the following are present
    in the same fused scaffold:
    1. Beta–lactam ring (4-membered cyclic amide): SMARTS "N1C(=O)CC1"
    2. Dimethyl substituents on an S–bound carbon (indicative of 2,2-dimethyl substitution): SMARTS "[S]-C([CH3])([CH3])"
    3. A carboxylate substituent (C(=O)O or C(=O)[O-])
    4. A carboxamido group (NC(=O))
    
    For each required functionality, we ensure that it is “close” (<= 4 bonds) to at least one atom
    of the beta–lactam ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as penicillin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If the molecule contains multiple fragments (e.g. salts), keep only the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Define SMARTS patterns for required substructures
    
    # (1) Beta-lactam ring: a 4-membered cyclic amide.
    beta_lactam_smarts = "N1C(=O)CC1"
    beta_lactam_q = Chem.MolFromSmarts(beta_lactam_smarts)
    beta_matches = mol.GetSubstructMatches(beta_lactam_q)
    if not beta_matches:
        return False, "Missing beta-lactam (4-membered cyclic amide) ring"
    # For our purpose, use the first match as the anchor.
    beta_anchor = set(beta_matches[0])
    
    # Helper function: given a substructure query and a distance threshold, check that at least
    # one of its match atoms is within 'max_dist' bonds of any atom in the anchor set.
    def is_close(query, max_dist=4):
        query_mol = Chem.MolFromSmarts(query)
        matches = mol.GetSubstructMatches(query_mol)
        if not matches:
            return False, f"Missing required group: {query}"
        # For each match, check distance (using shortest paths) to the beta anchor.
        for match in matches:
            for a in match:
                for anchor_atom in beta_anchor:
                    # Compute shortest path length between atoms a and anchor_atom.
                    sp = Chem.GetShortestPath(mol, a, anchor_atom)
                    # Number of bonds = (len(sp)-1)
                    if len(sp) - 1 <= max_dist:
                        return True, ""
        return False, f"Required group {query} not connected closely to beta-lactam ring"
    
    # (2) Dimethyl substituent on S-bound carbon: "[S]-C([CH3])([CH3])"
    dimethyl_smarts = "[S]-C([CH3])([CH3])"
    close, msg = is_close(dimethyl_smarts)
    if not close:
        return False, "Missing or mis-connected S-bound carbon with two methyl substituents (expected at position 2) " + msg
    
    # (3) Carboxylate substituent: allow protonated or deprotonated.
    carboxylate_smarts1 = "C(=O)[O-]"
    carboxylate_smarts2 = "C(=O)O"
    close1, msg1 = is_close(carboxylate_smarts1)
    close2, msg2 = is_close(carboxylate_smarts2)
    if not (close1 or close2):
        return False, "Missing carboxylate substituent (C(=O)O or C(=O)[O-]) near beta-lactam ring"
    
    # (4) Carboxamido group: "NC(=O)"
    carboxamido_smarts = "NC(=O)"
    close, msg = is_close(carboxamido_smarts)
    if not close:
        return False, "Missing carboxamido group (NC(=O)) connected to the beta-lactam scaffold " + msg
    
    return True, "Matches penicillin criteria: Contains beta‐lactam, 2,2-dimethyl, carboxylate and carboxamido substituents"

# For testing, one might run:
# test_smiles = "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)CCC[C@H](N)C(O)=O)"
# print(is_penicillin(test_smiles))