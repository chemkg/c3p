"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent on the alpha carbon.
Improved criteria:
  1. Identify a carboxylic acid group in a tolerant way (using SMARTS [CX3](=O)[OX2H,OX2-]).
  2. Ensure exactly one acid group is present.
  3. From the acid carbon, get candidate alpha carbons (neighbors via a single bond).
  4. For each candidate, accept it if:
       a. Its heavy-atom connectivity is “reasonable” (1–3 heavy neighbors).
       b. It has at least one double bond to an oxygen where that oxygen is “terminal”
          (i.e. oxygen’s degree is 1) and not part of the acid group.
If at least one candidate meets these criteria, the molecule is classified as a 2-oxo monocarboxylic acid.
Note: This is an ad hoc method and will not capture every nuance.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    
    Criteria:
      1. The molecule must contain exactly one carboxylic acid group (SMARTS: [CX3](=O)[OX2H,OX2-]).
      2. The carboxyl carbon must be attached (via a single bond) to at least one carbon (alpha carbon).
      3. At least one candidate alpha carbon must bear a ketone substituent (a double bond to oxygen)
         such that the oxygen is terminal (degree 1) and not also attached to the acid carbon.
      4. The alpha carbon is required to have a “reasonable” heavy-atom connectivity (1 to 3 heavy neighbors).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if it is classified as a 2-oxo monocarboxylic acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for the carboxylic acid group. We allow both protonated and deprotonated forms.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H,OX2-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    if not acid_matches:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, "More than one carboxylic acid group found; not a monocarboxylic acid"
    
    # Our SMARTS pattern gives us the acid carbon at index 0 in the match.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Identify candidate alpha carbons: they must be carbon atoms (atomic number = 6)
    # that are attached by a SINGLE bond to the acid carbon.
    candidate_alpha_atoms = []
    for nbr in acid_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            bond = mol.GetBondBetweenAtoms(acid_carbon_idx, nbr.GetIdx())
            if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                candidate_alpha_atoms.append(nbr)
    
    if not candidate_alpha_atoms:
        return False, "No alpha carbon (carbon neighbor) found attached to the carboxyl group"
    
    # For each candidate alpha carbon, check its connectivity and search for a ketone substituent.
    for alpha in candidate_alpha_atoms:
        # Count heavy-atom neighbors (exclude hydrogens)
        heavy_neighbors = [n for n in alpha.GetNeighbors() if n.GetAtomicNum() > 1]
        if not (1 <= len(heavy_neighbors) <= 3):
            # Skip candidates that are either under- or over-connected.
            continue
        
        # Look for a double bond to oxygen (C=O) that is not part of the acid group.
        for bond in alpha.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(alpha)
                if other.GetAtomicNum() == 8:  # candidate double-bonded oxygen
                    # Check that this oxygen is not part of the acid group.
                    # (i.e. it should not be connected to the acid carbon)
                    if acid_carbon_idx in [n.GetIdx() for n in other.GetNeighbors()]:
                        continue
                    # Additionally, require that the oxygen seems “terminal” (only one connection)
                    # which is typical for a ketone carbonyl oxygen.
                    if other.GetDegree() != 1:
                        continue
                    # If we have passed these tests, we have found the required 2-oxo motif.
                    return True, "Found a monocarboxylic acid with a 2-oxo substituent on the alpha carbon"
    
    # No candidate alpha carbon meets the criteria.
    return False, "No suitable alpha carbon with a terminal C=O substituent (2-oxo group) found on the acid group"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "C[C@H](C(=O)C(O)=O)c1c[nH]c2ccccc12",              # (S)-3-(indol-3-yl)-2-oxobutyric acid
        "Nc1ccc(O)c(\\C=C\\C(=O)C(O)=O)c1",                  # (3E)-4-(5-amino-2-hydroxyphenyl)-2-oxobut-3-enoic acid
        "O(C1=CC=C(CC(=O)C(O)=O)C=C1)C",                      # 3-(4-methoxyphenyl)-2-oxopropanoic acid
        "CC(C)(C)C(=O)C(O)=O",                               # 3,3-dimethyl-2-oxobutanoic acid
        "OC(=O)C(CCCCN)=O",                                  # 6-amino-2-oxohexanoic acid
        "CC(=O)C(O)=O",                                     # Pyruvic acid
        "CCCCCC(=O)C(O)=O",                                 # 2-oxooctanoic acid
        # False positives (should not classify):
        "F[C@H]1C[C@]2([C@]3([C@@]([C@H]([C@@H](C3)C)C(=O)C(O)=O)(C[C@H](O)[C@@]2([C@@]4(C1=CC(=O)C=C4)C)[H])C)[H])[H]",  # Fluocortin
        "OCC(=O)[C@@H](O)[C@H](O)C(=O)C(O)=O",             # 2,5-didehydro-D-gluconic acid
    ]
    
    for sm in test_smiles:
        result, reason = is_2_oxo_monocarboxylic_acid(sm)
        print(f"SMILES: {sm}")
        print(f"Result: {result}")
        print(f"Reason: {reason}")
        print("-" * 40)