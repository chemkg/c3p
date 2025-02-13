"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more free amino groups.
A free amino group is defined here as a sp3 nitrogen atom that has at least one hydrogen attached
and is not directly bonded to a carbonyl carbon (to avoid misclassifying amides).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is any organic compound that contains at least two free amino groups.
    A "free amino group" is defined here as a sp3 nitrogen (with at least one attached hydrogen)
    that is not directly bonded to a carbonyl group (i.e. not part of an amide bond).
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polyamine, False otherwise.
        str: Reason for the classification result.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Require an organic molecule: must contain at least one carbon atom
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain any carbon atoms, thus not organic"
    
    # Optional: if the molecule has amide bonds and is heavy, it might be a peptide/polyamide.
    # The simple SMARTS "C(=O)N" is used here.
    amide_query = Chem.MolFromSmarts("C(=O)N")
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol.HasSubstructMatch(amide_query) and mol_wt > 300:
        return False, "Contains amide bonds with high molecular weight, likely a peptide or polyamide"
    
    # Use a SMARTS query to count free amino groups.
    # The SMARTS "[NX3;H;!$([NX3][C](=O))]" means:
    #   [NX3]   -> a trivalent (sp3) nitrogen,
    #   H       -> that has at least one attached hydrogen,
    #   !$([NX3][C](=O)) -> and is not directly bonded to a carbon that is double-bonded to oxygen.
    free_amino_smarts = Chem.MolFromSmarts("[NX3;H;!$([NX3][C](=O))]")
    matches = mol.GetSubstructMatches(free_amino_smarts)
    # Extract the unique nitrogen atom indices from matches
    free_amino_indices = set(match[0] for match in matches)
    free_amino_count = len(free_amino_indices)
    
    if free_amino_count < 2:
        return False, f"Contains {free_amino_count} free amino group(s), need at least 2 for polyamine"
    else:
        return True, f"Contains {free_amino_count} free amino groups, satisfying polyamine criteria"

# Example usage (can be removed or commented out in production)
if __name__ == "__main__":
    # Test using trimethylenediamine as a simple example
    test_smiles = "NCCCN"  # trimethylenediamine
    result, reason = is_polyamine(test_smiles)
    print(result, reason)