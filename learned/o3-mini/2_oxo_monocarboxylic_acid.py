"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent on the alpha carbon.
Improved criteria:
  1. Identify a carboxylic acid group in a tolerant way (allowing [OX2H] or [OX2-]).
  2. Require that exactly one acid group is present.
  3. From the acid carbon, get candidate alpha carbons (neighbors that are carbon and share a single bond).
  4. For each candidate, allow a flexible heavy-atom connectivity (2–3 heavy neighbors) and then 
     check for at least one double bond to oxygen that is not part of the acid group.
     
Note: Because a number of these molecules exist in tautomeric or conjugated forms, this
processor uses an ad hoc method that may not capture every nuance. If no candidate passes the test,
the molecule is not classified.
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    
    Criteria:
      1. The molecule must contain exactly one carboxylic acid group defined as [CX3](=O)[OX2H,OX2-].
      2. The carboxyl carbon must be attached (via a single bond) to at least one carbon (the alpha carbon).
      3. The candidate alpha carbon must have a double bond to an oxygen atom (the keto group)
         that is not simply part of the acid group.
      4. To reduce mis‐classifications in complex molecules, we require that the candidate alpha
         carbon has a “reasonable” heavy‐atom connectivity (2 or 3 neighbors).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a 2-oxo monocarboxylic acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define the SMARTS pattern for a carboxylic acid group.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H,OX2-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    if not acid_matches:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, "More than one carboxylic acid group found; not a monocarboxylic acid"
        
    # Get the acid group. By our SMARTS, atom index 0 is the acid (carboxyl) carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Look for candidate alpha carbons: they must be carbons directly attached by a SINGLE bond
    candidate_alpha_atoms = []
    for nbr in acid_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # candidate must be carbon
            bond = mol.GetBondBetweenAtoms(acid_carbon_idx, nbr.GetIdx())
            if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                candidate_alpha_atoms.append(nbr)
    
    if not candidate_alpha_atoms:
        return False, "No alpha carbon (carbon neighbor) found attached to the carboxyl group"
        
    # For each candidate alpha carbon, check connectivity and for a keto (C=O) substituent.
    for alpha_atom in candidate_alpha_atoms:
        # Allow a flexible heavy-atom neighbor count: typically expect 2 or 3 heavy neighbors.
        heavy_neighbors = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) < 2 or len(heavy_neighbors) > 3:
            # Skip over-elaborated or under-connected candidates.
            continue
        
        # Check each bond from the alpha carbon for a double bond to an oxygen.
        for bond in alpha_atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other_atom = bond.GetOtherAtom(alpha_atom)
                if other_atom.GetAtomicNum() == 8:
                    # If this oxygen is bonded to the acid carbon (i.e. part of the acid group) then skip.
                    if acid_carbon_idx in [nbr.GetIdx() for nbr in other_atom.GetNeighbors()]:
                        continue
                    # We have found a candidate: the alpha carbon bears a C=O group.
                    return True, "Found a monocarboxylic acid with a 2-oxo substituent on the alpha carbon"
    
    return False, "No suitable alpha carbon with a C=O substituent (2-oxo group) found on the acid group"

# Example usage (testing a few representative cases):
if __name__ == "__main__":
    test_smiles = [
        # True positives
        "C[C@H](C(=O)C(O)=O)c1c[nH]c2ccccc12",              # (S)-3-(indol-3-yl)-2-oxobutyric acid
        "Nc1ccc(O)c(\\C=C\\C(=O)C(O)=O)c1",                  # (3E)-4-(5-amino-2-hydroxyphenyl)-2-oxobut-3-enoic acid
        "O(C1=CC=C(CC(=O)C(O)=O)C=C1)C",                      # 3-(4-methoxyphenyl)-2-oxopropanoic acid
        "CC(C)(C)C(=O)C(O)=O",                               # 3,3-dimethyl-2-oxobutanoic acid
        "OC(=O)C(CCCCN)=O",                                  # 6-amino-2-oxohexanoic acid
        "CC(=O)C(O)=O",                                     # pyruvic acid
        # Some examples that were previously false positives or negatives:
        "F[C@H]1C[C@]2([C@]3([C@@]([C@H]([C@@H](C3)C)C(=O)C(O)=O)(C[C@H](O)[C@@]2([C@@]4(C1=CC(=O)C=C4)C)[H])C)[H])[H]",  # Fluocortin (should not classify)
        "O(C=1C=C(C=CC1)/C=C(\\O)/C(O)=O)C",                # (2Z)-2-hydroxy-3-(4-methoxyphenyl)prop-2-enoic acid (challenging case)
    ]
    
    for sm in test_smiles:
        result, reason = is_2_oxo_monocarboxylic_acid(sm)
        print(f"SMILES: {sm}")
        print(f"Result: {result}")
        print(f"Reason: {reason}")
        print("-" * 40)