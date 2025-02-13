"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: Tricarboxylic acid – An oxoacid containing three carboxy groups.
This classifier not only counts the free (acidic) carboxyl groups but also
applies extra filters. The idea is that a “true” tricarboxylic acid (such as citric acid,
nitrilotriacetic acid, aconitic acid, etc.) usually is not very complex and does
not contain multiple amide bonds (which are present in peptides or conjugated metabolites).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    
    A tricarboxylic acid is defined here as a molecule with exactly three “free”
    carboxyl groups (i.e. in either the protonated form C(=O)[OH] or the deprotonated form C(=O)[O-])
    that is not otherwise excessively complex (for example containing multiple amide bonds
    or having a large heavy-atom count) as occurs in peptides or large glycoconjugates.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True if molecule is classified as a tricarboxylic acid with a reason;
                   False otherwise, with the reason for rejection.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for a carboxyl group in protonated (acid) and deprotonated forms.
    # The patterns match the carbonyl carbon bonded to an -OH or -O(-) group.
    carboxyl_neutral = Chem.MolFromSmarts("C(=O)[OX2H1]")
    carboxyl_anion   = Chem.MolFromSmarts("C(=O)[OX1-]")
    
    # Get all unique matches.
    matches_neutral = set(mol.GetSubstructMatches(carboxyl_neutral))
    matches_anion   = set(mol.GetSubstructMatches(carboxyl_anion))
    carboxyl_matches = matches_neutral.union(matches_anion)
    n_carboxyl = len(carboxyl_matches)
    
    if n_carboxyl != 3:
        return False, f"Found {n_carboxyl} carboxyl group(s), but exactly 3 are required."
    
    # Compute some overall properties.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    num_rings = mol.GetRingInfo().NumRings()
    
    # Heuristic 1: Exclude overly complex molecules.
    # True tricarboxylic acids tend to be relatively “small/simple.”
    # (Thresholds here are heuristic; many genuine oxoacids have heavy atom counts below 60 and few rings.)
    if heavy_atoms > 60:
        return False, "Molecule is very complex (more than 60 heavy atoms) and unlikely to be a simple oxoacid."
    if num_rings > 5:
        return False, "Molecule contains many rings and appears too complex for a typical tricarboxylic acid."
    
    # Heuristic 2: Exclude cases with multiple amide bonds.
    # We use a simple SMARTS for an amide bond.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    n_amide = len(amide_matches)
    # Allow one amide bond (which can occur in small molecules like nitrilotriacetic acid derivatives)
    # but if there are two or more and the molecular weight is relatively high, we assume a peptide or conjugate.
    if n_amide >= 2 and mol_wt > 250:
        return False, "Molecule contains multiple amide bonds, suggesting a peptide-like or conjugated structure."
    
    # Passed all tests.
    return True, "Molecule contains exactly 3 free carboxyl groups and meets additional structural filters for an oxoacid."

# Example usage (uncomment for testing):
# test_smiles = "OC(=O)CC(O)(CC(O)=O)C(O)=O"  # Example: citric acid
# result, reason = is_tricarboxylic_acid(test_smiles)
# print(result, reason)