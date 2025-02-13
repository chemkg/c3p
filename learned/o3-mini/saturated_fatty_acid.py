"""
Classifies: CHEBI:26607 saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Saturated fatty acid

A saturated fatty acid is defined as a fatty acid that contains exactly one carboxylic acid group,
has a long aliphatic chain containing no carbon–carbon multiple bonds (apart from the C=O of the acid),
and consists only of carbon, hydrogen, and oxygen atoms.
It is known to produce adverse biological effects when ingested to excess.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    The molecule must have:
      - Only C, H, and O atoms.
      - Exactly one carboxylic acid group (using the pattern [CX3](=O)[OX2H]).
      - No carbon-carbon double or triple bonds except for the carbonyl in the acid group.
      - No rings (to ensure the structure is an open-chain fatty acid).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a saturated fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule contains only allowed atoms: C (6), O (8), and H (1)
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Contains atoms other than C, H, and O"
    
    # Check for rings: a fatty acid should not contain any rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structures"
    
    # Define a SMARTS pattern for the carboxylic acid group (C(=O)[OH])
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    if len(acid_matches) != 1:
        return False, "More than one carboxylic acid group found"
    
    # Check for carbon-carbon multiple bonds (i.e. C=C or C≡C) in parts of the molecule 
    # other than the carbonyl group of the acid.
    # We iterate over bonds, and if a bond is UNSATURATED (DOUBLE or TRIPLE) and both atoms are carbon,
    # then it is not a saturated fatty acid.
    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        # If the bond is double or triple
        if bond_type in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # ignore the carbonyl double bond because one of the oxygens is involved
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                return False, "Contains carbon-to-carbon multiple bonds"
    
    return True, "Molecule has a single carboxylic acid group, only C/H/O atoms, no rings, and no C-C unsaturation"

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "OCCCCCCCCCCCCCCCCCCCC(O)=O",  # 20-hydroxyicosanoic acid (saturated fatty acid)
        "S1N(C2=CC=C(C=C2)C(OCC)=O)C(=O)C=C1"  # false example with non C/H/O atoms and aromatic rings
    ]
    for sm in test_smiles:
        result, reason = is_saturated_fatty_acid(sm)
        print(f"SMILES: {sm} -> {result}, Reason: {reason}")