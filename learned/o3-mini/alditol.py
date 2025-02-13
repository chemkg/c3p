"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol 
Definition: A carbohydrate that is an acyclic polyol having the general formula 
HOCH2[CH(OH)]nCH2OH (formally derivable from an aldose by reduction of the carbonyl group).

This version uses a more refined SMARTS query by explicitly requiring:
  - Terminal carbons: [CH2;!R]([OX2H])  (a non‐ring carbon with two hydrogens bonded to –OH)
  - Interior carbons: [CH;!R]([OX2H])  (a non‐ring carbon with one hydrogen bonded to –OH)
and connecting them via single bonds.
"""

from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule can be classified as an alditol 
    based on its SMILES string.
    
    The approach is as follows:
      1) Reject the molecule if it contains any carbonyl group ([#6]=[OX1]).
      2) For a series of chain lengths (n = number of carbons from 3 to 12),
         build a SMARTS pattern corresponding to:
           Terminal carbons: HOCH2 group -> "[CH2;!R]([OX2H])"
           Interior carbons: CH(OH) groups -> "[CH;!R]([OX2H])"
         Referring to the general pattern HOCH2[CH(OH)]^(n-2)CH2OH.
      3) Attempt a substructure match for these patterns.
      4) If any match is found, then we classify the structure as an alditol.
    
    Args:
        smiles (str): The SMILES string describing the molecule.
    
    Returns:
        bool: True if a matching acyclic polyol chain is found, False otherwise.
        str: A reason describing the result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens to standardize neighbor counts
    mol = Chem.AddHs(mol)
    
    # Rule out molecules that contain a carbonyl group,
    # since alditols must be fully reduced (no C=O present).
    carbonyl_query = Chem.MolFromSmarts("[#6]=[OX1]")
    if mol.HasSubstructMatch(carbonyl_query):
        return False, "Contains a carbonyl group, so not a fully reduced (alditol) structure"
    
    # Define SMARTS for the terminal and interior units:
    # Terminal carbon (CH2OH): non-ring CH2 bonded to an -OH.
    term_smarts = "[CH2;!R]([OX2H])"
    # Interior carbon (CH(OH)): non-ring CH bonded to an -OH.
    interior_smarts = "[CH;!R]([OX2H])"
    
    # Try chain lengths from 3 to 12 carbons.
    for n in range(3, 13):
        # Build SMARTS string for a linear chain:
        # Example for n=3: terminal-interior-terminal
        pattern_parts = []
        # first terminal:
        pattern_parts.append(term_smarts)
        # (n-2) interior units
        for _ in range(n - 2):
            pattern_parts.append(interior_smarts)
        # last terminal:
        pattern_parts.append(term_smarts)
        # Join with a single-bond connector (the '-' indicates a single bond)
        pattern_smarts = "-".join(pattern_parts)
        
        query = Chem.MolFromSmarts(pattern_smarts)
        if query is None:
            continue  # should not happen, but safety first
        
        # Get substructure matches; not using chirality constraints for flexibility.
        matches = mol.GetSubstructMatches(query, useChirality=False)
        if matches:
            return True, f"Contains an acyclic polyol chain matching HOCH2[CH(OH)]^{n-2}CH2OH pattern"
    
    return False, "No acyclic polyol chain matching the required HOCH2[CH(OH)]nCH2OH pattern was found"

# If run as a script, do simple tests.
if __name__ == "__main__":
    # Example test SMILES for alditols.
    test_smiles = [
        "OC[C@H](O)[C@H](O)CO",        # erythritol
        "OC(C(O)C(O)C(O)CO)C(O)C(O)CO", # D-Erythro-D-galacto-octitol
        "OCC(O)CO",                    # glycerol (although a polyol, it may be too short)
        "OC[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO"  # D-iditol
    ]
    for smi in test_smiles:
        result, reason = is_alditol(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")