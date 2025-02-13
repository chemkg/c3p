"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
#!/usr/bin/env python
"""
Classifies: Tricarboxylic acid – An oxoacid containing three carboxy groups.
This classifier identifies “free” carboxyl groups by matching SMARTS for the neutral 
and deprotonated forms and then confirming that each carboxyl group is terminal in the sense
that its carbon is only connected to exactly one non-oxygen atom (and no nitrogen, which might indicate 
an amide linkage). We then check that the molecule contains exactly three such groups and that the 
density of these groups relative to the heavy-atom count is high enough (to avoid cases where 3 carboxyl groups 
are “diluted” in an otherwise complex molecule).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    
    A tricarboxylic acid is defined here as an oxoacid containing exactly three 'free' carboxyl groups.
    In our approach, a "free" carboxyl group is detected via SMARTS (for both protonated and deprotonated forms)
    and further validated by checking that the carboxyl carbon is bonded to exactly two oxygens (the carbonyl 
    and hydroxyl/oxide) and one other atom that is not oxygen or nitrogen. This latter check aims to rule out carboxyl
    groups involved in amide bonds.
    
    Furthermore, a relative “density” filter is applied: if the overall heavy atom count is high, the proportion
    of free carboxyl groups must also be high – a property expected for genuine simple oxoacids but not for larger, 
    conjugated molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a tricarboxylic acid, False otherwise.
        str: An explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for a free carboxyl group.
    # We use explicit atom maps to ensure that the first atom is the carboxyl carbon.
    carboxyl_neutral = Chem.MolFromSmarts("[C:1](=O)[O;X2H1]")
    carboxyl_anion   = Chem.MolFromSmarts("[C:1](=O)[O;X1-]")
    
    # Get matches for both protonation states.
    matches_neutral = mol.GetSubstructMatches(carboxyl_neutral)
    matches_anion   = mol.GetSubstructMatches(carboxyl_anion)
    
    # Use a set to keep track of unique carboxyl group carbon atom indices.
    free_carboxyl_c_indices = set()
    
    # Define a helper function to check if a matched carboxyl group is “free.”
    # A free carboxyl group, for our purposes, is one where the carboxyl carbon (match index 0)
    # has exactly three bonds: two to oxygen and one to some other atom (the R group).
    # Furthermore, if that non-oxygen neighbor is nitrogen, we suspect it may be part of an amide.
    def is_free_carboxyl(match):
        c_idx = match[0]
        atom_c = mol.GetAtomWithIdx(c_idx)
        neighbors = atom_c.GetNeighbors()
        # Count non-oxygen neighbors.
        non_o_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() != 8]
        # If any neighbor is nitrogen (atomic number 7), consider it not free.
        if any(nbr.GetAtomicNum() == 7 for nbr in non_o_neighbors):
            return False
        # For a typical carboxyl carbon, it should have exactly one non-oxygen neighbor.
        if len(non_o_neighbors) != 1:
            return False
        # Also ensure that the total degree is exactly 3 (two oxygens and one R-group).
        if atom_c.GetDegree() != 3:
            return False
        return True
    
    # Check each match in the neutral and anionic forms.
    for match in matches_neutral:
        if is_free_carboxyl(match):
            free_carboxyl_c_indices.add(match[0])
    for match in matches_anion:
        if is_free_carboxyl(match):
            free_carboxyl_c_indices.add(match[0])
    
    n_free_carboxyl = len(free_carboxyl_c_indices)
    
    # We require exactly three free carboxyl groups.
    if n_free_carboxyl != 3:
        return False, f"Found {n_free_carboxyl} free carboxyl group(s); exactly 3 are required."
    
    # Compute overall heavy atom count.
    num_heavy = mol.GetNumHeavyAtoms()
    
    # Apply a density heuristic: in genuine (often small/simple) tricarboxylic acids
    # the three carboxyl groups represent a significant fraction of the heavy atoms.
    density = n_free_carboxyl / num_heavy  # ratio of free carboxyl groups to heavy atoms
    # In our tests, genuine tricarboxylic acids typically have high ratios (e.g., above 0.1).
    # If the molecule is quite heavy and the free acid groups are diluted (density < 0.1),
    # we suspect the overall structure is not a simple oxoacid.
    if num_heavy >= 40 and density < 0.1:
        return False, ("Molecule has 3 free carboxyl groups but a low carboxyl-group density "
                       f"({density:.2f}) given its {num_heavy} heavy atoms, suggesting it is an overly complex structure.")
    
    # Passed all tests.
    return True, ("Molecule contains exactly 3 free carboxyl groups with sufficient density "
                  f"(density = {density:.2f}, heavy atoms = {num_heavy}), consistent with a tricarboxylic acid.")

# Example usage (for testing, uncomment the following lines):
# test_smiles = "OC(=O)CC(O)(CC(O)=O)C(O)=O"  # Example: citric acid
# result, reason = is_tricarboxylic_acid(test_smiles)
# print(result, reason)