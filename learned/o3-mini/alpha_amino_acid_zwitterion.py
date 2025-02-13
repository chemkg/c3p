"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid-zwitterion obtained by transfer of a proton from the carboxy to the amino group 
of any alpha-amino acid; major species at pH 7.3.
The key motif is an alpha-carbon (with at least one hydrogen) directly bonded to a positively charged nitrogen 
and a carboxylate group.
"""

from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    The method looks for an alpha-carbon that is bound to a positively charged amino group and a carboxylate.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule contains an alpha-amino-acid zwitterion core, otherwise False.
        str: A message explaining the basis for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS for the alpha-amino-acid zwitterion core.
    # This pattern matches a tetrahedral (sp3) carbon with one explicit hydrogen ([C;H1]) 
    # that is bonded to a positively charged nitrogen ([N+]) and a carboxylate moiety (C(=O)[O-]).
    # The order of substituents in the SMARTS does not force the molecule to be drawn in that order,
    # but rather ensures that the three groups are attached to the candidate alpha carbon.
    aa_zwitterion_pattern = Chem.MolFromSmarts("[C;H1]([N+])C(=O)[O-]")
    
    # Check for the presence of the alpha amino acid zwitterion substructure.
    if mol.HasSubstructMatch(aa_zwitterion_pattern):
        return True, "Molecule contains an alpha-amino-acid zwitterion core"
    else:
        return False, "Alpha-amino-acid zwitterion core not found"

# (Optional) if you wish to test with one of the provided examples, you can run:
# test_smiles = "CC(C)[C@H]([NH3+])C([O-])=O"  # L-valine zwitterion
# print(is_alpha_amino_acid_zwitterion(test_smiles))