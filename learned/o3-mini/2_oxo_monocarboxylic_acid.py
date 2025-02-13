"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent on the alpha carbon.
Improved criteria:
  1. Identify a carboxylic acid group in a tolerant way (allow [OX2H] or [OX2-]).
  2. Require that exactly one acid group is present.
  3. From the acid carbon, get non‐oxygen neighbors (alpha carbons). For each candidate,
     check that it has a double bond to an oxygen (that is not the carboxyl oxygen)
     and that its heavy atom connectivity is consistent with a simple 2‑oxo motif (exactly 3 neighbors).
     
Note: This ad hoc method may not capture every nuance of challenging natural product structures,
but it attempts to reduce both false positives and false negatives.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    The criteria are:
      1. The molecule must have exactly one carboxylic acid group.
         (The acid group is identified by [CX3](=O)[OX2H,OX2-]).
      2. The carboxyl carbon must be attached via a single bond to an alpha carbon.
      3. That alpha carbon (i.e. the one directly bonded to the acid carbon) must have a 
         keto (C=O) group attached. In order to reduce mis‐classifications, we also check that
         the alpha carbon is not “over‐substituted” (i.e. it exactly connects to three heavy atoms:
         the acid carbon, the keto oxygen and one other substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a 2-oxo monocarboxylic acid, False otherwise.
        str: Explanation of the reasoning.
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carboxylic acid group.
    # This pattern allows for either a hydroxyl (-OH) or a deprotonated acid (-O^-)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H,OX2-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    if len(acid_matches) == 0:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, "More than one carboxylic acid group found; not a monocarboxylic acid"
    
    # Get the acid match: by our SMARTS the first atom (index 0) is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Find the candidate alpha carbon as any neighbor that is not oxygen.
    candidate_alpha_atoms = []
    for neighbor in acid_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() != 8:  # skip oxygens (part of the acid group)
            candidate_alpha_atoms.append(neighbor)
    
    if not candidate_alpha_atoms:
        return False, "No alpha carbon (non-oxygen neighbor) found attached to the carboxyl group"
    
    # For each candidate alpha atom, attempt to confirm it carries the 2-oxo (ketone) substituent.
    for alpha_atom in candidate_alpha_atoms:
        # Check heavy atom connectivity of alpha: expect exactly 3 neighbors
        heavy_neighbors = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # Our expectation: one bond is to the acid carbon, one should be the keto oxygen and one other group.
        if len(heavy_neighbors) != 3:
            # if not exactly three, skip this candidate (might be part of complex ring systems)
            continue

        # Now check that one of the bonds from alpha_atom is a C=O bond (not the acid bond).
        has_oxo = False
        for bond in alpha_atom.GetBonds():
            # We want a double bond to oxygen.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other_atom = bond.GetOtherAtom(alpha_atom)
                # skip if this oxygen is already connected to the acid carbon 
                # (i.e. one of the acid group's oxygens)
                if other_atom.GetAtomicNum() == 8:
                    # We double-check: if this oxygen is also bonded to the acid carbon then it might be part of the acid group.
                    # So we exclude it if its neighbor list includes our acid carbon.
                    if acid_carbon_idx in [nbr.GetIdx() for nbr in other_atom.GetNeighbors()]:
                        continue
                    # Otherwise, this double bond qualifies as the 2-oxo substituent.
                    has_oxo = True
                    break
        
        if has_oxo:
            return True, "Found a monocarboxylic acid with a 2-oxo substituent on the alpha carbon"
    
    # If no candidate alpha atom meets the criteria, classification fails.
    return False, "No suitable alpha carbon with a C=O substituent (2-oxo group) found on the acid group"

# Example usage:
if __name__ == "__main__":
    # List of example SMILES (both from true positives and edge cases)
    examples = [
        "C[C@H](C(=O)C(O)=O)c1c[nH]c2ccccc12",  # true positive
        "Nc1ccc(O)c(\\C=C\\C(=O)C(O)=O)c1",         # true positive
        "F[C@H]1C[C@]2([C@]3([C@@]([C@H]([C@@H](C3)C)C(=O)C(O)=O)(C[C@H](O)[C@@]2([C@@]4(C1=CC(=O)C=C4)C)[H])C)[H])[H]",  # false positive: fluocortin
        "OC(=O)C(CCCCN)=O",                         # true positive: 6-amino-2-oxohexanoic acid
        "CC(=O)C(O)=O",                            # true positive: pyruvic acid
    ]
    
    for sm in examples:
        result, reason = is_2_oxo_monocarboxylic_acid(sm)
        print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n{'-'*40}")