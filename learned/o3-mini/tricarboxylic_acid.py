"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
#!/usr/bin/env python
"""
Classifies: Tricarboxylic acid – An oxoacid containing three carboxy groups.
This classifier identifies “free” carboxyl groups by matching SMARTS for the neutral 
and deprotonated forms and then confirming that the carboxyl carbon is terminal 
(i.e. it is attached to exactly two oxygens and one other atom that is not oxygen or nitrogen).
Finally, the code computes the “density” as the number of free carboxyl groups divided by 
the molecule’s heavy atom count. Because genuine (often small/simple) oxoacids have the acid groups 
as a significant fraction of their heavy atoms, we use a piecewise threshold:
    • For molecules with <20 heavy atoms, require density >= 0.20.
    • For molecules with 20–40 heavy atoms, require density >= 0.11.
    • For molecules with ≥40 heavy atoms, require density >= 0.10.
These heuristics help avoid classifying peptide-like (or overly complex) molecules as tricarboxylic acids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid (an oxoacid with exactly 3 “free” carboxyl groups)
    based on its SMILES string.
    
    A free carboxyl group is detected in both its protonated [C](=O)[O;H1] and deprotonated forms [C](=O)[O-],
    and then further required to be terminal; namely, the carboxyl carbon must be attached to exactly 3 atoms,
    of which 2 must be oxygens (for the carbonyl and hydroxyl/oxide) and the one remaining, the R-group, must 
    not be oxygen or nitrogen (to avoid, e.g., amide bonds). Finally, a density check is applied: the ratio 
    (free carboxyl groups)/(heavy atoms) must be above a threshold that depends on the total heavy atom count.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a tricarboxylic acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for free carboxyl groups
    # The pattern ensures that we are matching a carboxyl carbon attached (via double bond) to oxygen,
    # and to an –OH (or –O– in the anion form). We use an atom map to later identify the carboxyl carbon.
    carboxyl_neutral = Chem.MolFromSmarts("[C:1](=O)[O;X2H1]")
    carboxyl_anion   = Chem.MolFromSmarts("[C:1](=O)[O;X1-]")

    # Find matches for both forms.
    matches_neutral = mol.GetSubstructMatches(carboxyl_neutral)
    matches_anion   = mol.GetSubstructMatches(carboxyl_anion)

    # Use a set to store the unique index of the carboxyl carbon for each free carboxyl group.
    free_carboxyl_c_indices = set()

    # Helper function to decide if a matched carboxyl group is “free.”
    # Conditions:
    #   - The carboxyl carbon (first atom in the match) must have exactly 3 bonds.
    #   - Two of its bonds must be to oxygen.
    #   - Its one non-oxygen neighbor (the R group) should not be nitrogen.
    def is_free_carboxyl(match):
        c_idx = match[0]
        atom_c = mol.GetAtomWithIdx(c_idx)
        # Check total degree: for a carboxyl carbon, expect 3 connections (C=O, C-O and R group)
        if atom_c.GetDegree() != 3:
            return False
        neighbors = atom_c.GetNeighbors()
        # Count oxygen neighbors
        oxygen_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
        # There must be exactly 2 oxygens (one double-bonded carbonyl and one –OH or –O–)
        if len(oxygen_neighbors) != 2:
            return False
        # Identify the one non-oxygen neighbor (the R-group)
        non_o_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() != 8]
        if len(non_o_neighbors) != 1:
            return False
        # If the R-group is nitrogen, then it might be part of an amide or related motif
        if non_o_neighbors[0].GetAtomicNum() == 7:
            return False
        return True

    # Evaluate matches from both SMARTS patterns.
    for match in matches_neutral:
        if is_free_carboxyl(match):
            free_carboxyl_c_indices.add(match[0])
    for match in matches_anion:
        if is_free_carboxyl(match):
            free_carboxyl_c_indices.add(match[0])

    n_free = len(free_carboxyl_c_indices)

    # Require exactly 3 free carboxyl groups.
    if n_free != 3:
        return False, f"Found {n_free} free carboxyl group(s); exactly 3 are required."

    # Calculate the number of heavy atoms.
    num_heavy = mol.GetNumHeavyAtoms()
    # Calculate carboxyl group density.
    density = n_free / num_heavy

    # Apply a piecewise density threshold:
    #   For very small molecules (<20 heavy atoms): require density >= 0.20.
    #   For mid-sized molecules (20 to <40 heavy atoms): require density >= 0.11.
    #   For larger molecules (>=40 heavy atoms): require density >= 0.10.
    if num_heavy < 20:
        density_threshold = 0.20
    elif num_heavy < 40:
        density_threshold = 0.11
    else:
        density_threshold = 0.10

    if density < density_threshold:
        return False, (f"Molecule has 3 free carboxyl groups but a low carboxyl-group density "
                       f"({density:.2f}) given its {num_heavy} heavy atoms; expected at least {density_threshold:.2f} "
                       "for a genuine tricarboxylic acid (suggesting an overly complex structure).")

    # Passed all tests.
    return True, (f"Molecule contains exactly 3 free carboxyl groups with sufficient density "
                  f"(density = {density:.2f}, heavy atoms = {num_heavy}), consistent with a tricarboxylic acid.")

# Example usage:
# test_smiles = "OC(=O)CC(O)(CC(O)=O)C(O)=O"  # Example: citric acid
# result, explanation = is_tricarboxylic_acid(test_smiles)
# print(result, explanation)