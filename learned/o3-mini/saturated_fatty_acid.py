"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: Saturated fatty acid
A saturated fatty acid is defined as a fatty acid (i.e. contains a carboxylic acid group)
that does not contain any carbon-carbon multiple bonds.
Known to produce adverse biological effects when ingested to excess.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    A saturated fatty acid must have:
      - A carboxylic acid group (i.e. C(=O)O substructure).
      - No carbon to carbon double or triple bonds (unsaturation) in the hydrocarbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a saturated fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for the presence of a carboxylic acid group.
    # Using a simple SMARTS for carboxyl: this matches the pattern C(=O)O.
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Iterate over all bonds to detect any carbon-carbon multiple bonds.
    # Allow multiple bonds only if one of the atoms is not carbon (e.g., the C=O in COOH).
    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Check if both atoms are carbons
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                return False, "Contains carbon-to-carbon multiple bonds"
    
    return True, "Molecule has a carboxylic acid group and no carbon-carbon unsaturation"