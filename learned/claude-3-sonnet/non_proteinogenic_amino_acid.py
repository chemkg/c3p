"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI:33669 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is any amino acid that is not naturally encoded in the genetic code of any organism.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for amino acid backbone (-NH-CH(R)-COOH)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4;H2]([CH1]([CH3X4,CH2X3,CH1X2,CH0X1])C(=O)[OH,OX1])")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"
    
    # Check for proteinogenic amino acids
    proteinogenic_smarts = ["NC(C)C(=O)O", "NC(CC(=O)O)C(=O)O", "NC(CS)C(=O)O", "NC(Cc1ccccc1)C(=O)O",
                             "NC(Cc1c[nH]cn1)C(=O)O", "NC(CC(=O)O)C(=O)N", "NC(CC[C@H](C(=O)O)N)C(=O)O",
                             "NC(Cc1cnc[nH]1)C(=O)O", "NC(CC1=CNC2=C1C=CC=C2)C(=O)O", "NC(CC(N)=O)C(=O)O",
                             "NC(CC1=CC=CC=N1)C(=O)O", "NC(CC1=CN=CN1)C(=O)O", "NC(Cc1c[nH]c2c1cccc2)C(=O)O",
                             "NC(CO)C(=O)O", "NC(Cc1c[nH]c2ccccc12)C(=O)O", "NC(C(=O)O)C(C)(C)C",
                             "NC(Cc1ccc(O)cc1)C(=O)O", "NC(Cc1c[nH]cn1)C(=O)N", "NC(CCC(=O)O)C(=O)O"]
    for smarts in proteinogenic_smarts:
        proteinogenic_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(proteinogenic_pattern):
            return False, "Molecule is a proteinogenic amino acid"
    
    # Check for substitutions on the alpha carbon
    alpha_carbon_pattern = Chem.MolFromSmarts("[CH1]([CH3X4,CH2X3,CH1X2,CH0X1])[NX3,NX4;H2][CH3X4,CH2X3,CH1X2,CH0X1]")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if len(alpha_carbon_matches) == 0:
        return False, "No alpha carbon found in amino acid backbone"
    alpha_carbon_idx = alpha_carbon_matches[0][0]
    alpha_carbon_atom = mol.GetAtomWithIdx(alpha_carbon_idx)
    if sum(1 for a in alpha_carbon_atom.GetNeighbors() if a.GetAtomicNum() != 1) > 3:
        return True, "Contains a substituent on the alpha carbon"
    
    # Check for ring systems or cyclic structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return True, "Contains a ring system or cyclic structure"
    
    # Check for unusual functional groups or substitutions
    unusual_smarts = ["[N+](=O)[O-]", "[N+](=O)[O-]c", "[P+]", "Br", "Cl", "I", "F", "[O-]P", "[O-]P(=O)", "[O-]P(O)(O)=O"]
    for smarts in unusual_smarts:
        unusual_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(unusual_pattern):
            return True, f"Contains unusual functional group or substitution ({smarts})"
    
    # Count number of amino groups (-NH2) and carboxyl groups (-COOH)
    amino_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N'
                      and sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in atom.GetNeighbors()) == 1)
    carboxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C'
                         and atom.GetFormalCharge() == 0
                         and sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in atom.GetNeighbors()) == -1)
    if amino_count != 1 or carboxyl_count != 1:
        return False, "Must have exactly one amino and one carboxyl group"
    
    return True, "Contains an amino acid backbone with additional modifications or substituents"