"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: Alkanethiol
An alkanethiol is defined as a compound in which a sulfanyl group (-SH)
is attached to an alkyl group.
The function checks that:
  • the molecule has at least one –SH group (thiol) attached to an sp3
    (saturated) carbon,
  • the molecule does not have formal charges on non-counterion atoms,
  • the molecule does not contain amide bonds (suggestive of peptides/proteins),
  • the molecule is not overly complex (using a ring count/heavy atom count as proxy).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkanethiol(smiles: str):
    """
    Determine if a molecule (given by its SMILES) is classified as an alkanethiol.
    An alkanethiol is defined as having one or more -SH groups attached to sp3 (alkyl) carbons.
    
    Steps:
      1. Parse the SMILES and add explicit hydrogens.
      2. Exclude molecules with nonzero formal charges (except for typical counterions).
      3. Find all –SH groups using the SMARTS "[#16;X2&H1]-[C]".
         For each match, ensure that the carbon is sp3 hybridized.
      4. Exclude molecules that contain amide bonds (SMARTS "[NX3][CX3](=O)").
      5. Use ring count and heavy atom count as a proxy for over-complexity.
    
    Args:
        smiles (str): A SMILES string.
        
    Returns:
        bool: True if molecule is classified as alkanethiol; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the -SH group is properly represented.
    mol = Chem.AddHs(mol)
    
    # Exclude molecules carrying nonzero formal charges (ignoring common counterions)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ["Cl", "Br", "F", "I"]:
            if atom.GetFormalCharge() != 0:
                return False, "Molecule has non-zero formal charges (excluding typical counterions)"
    
    # Define the SMARTS for a thiol (-SH) attached to a carbon.
    thiol_smarts = "[#16;X2&H1]-[C]"
    thiol_pattern = Chem.MolFromSmarts(thiol_smarts)
    if thiol_pattern is None:
        return False, "Invalid thiol SMARTS pattern"
    
    # Find all occurrences of the thiol pattern.
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)
    if not thiol_matches:
        return False, "No thiol (-SH) group attached to a carbon found"
    
    # Check if at least one found thiol has the attached carbon in an sp3 environment.
    valid_thiol_found = False
    for match in thiol_matches:
        # Each match is a tuple (S_index, C_index)
        s_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        if c_atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            valid_thiol_found = True
            break
    if not valid_thiol_found:
        return False, "Thiol group found but not attached to an sp3 (alkyl) carbon"
    
    # Exclude molecules containing amide bonds (common in peptides/proteins).
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Molecule contains amide bonds, suggesting a peptide or protein derivative"
    
    # Use ring and heavy atom count to filter overly complex molecules.
    # Using rdMolDescriptors.CalcNumRings which returns an int.
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if ring_count > 1 and heavy_atom_count > 20:
        return False, "Molecule is too complex (multiple rings and many heavy atoms) to be a simple alkanethiol"
    
    # If we reached here, the molecule is a valid alkanethiol.
    return True, "Contains at least one -SH group attached to an sp3 alkyl carbon in a simple molecular context"

# Example test cases to verify functionality.
if __name__ == "__main__":
    test_cases = [
        ("Cl.SCCN", "Cysteamine hydrochloride"),
        ("CCS", "(ethylsulfanyl)methanethiol"),
        ("SCCCCCCS", "1,6-Hexanedithiol"),
        ("SC(C(C)C=O)CC", "3-Mercapto-2-methylpentanal"),
        ("SCCCCCC", "1-Hexanethiol"),
        ("S\\C=C(\\CC)/C", "2-Methyl-1-propenethiol"),
        ("C[C@H](O)[C@@H](S)P(O)(O)=O", "[(1R,2S)-2-hydroxy-1-thiopropyl]phosphonic acid"),
        ("CCSCCS", "2-(ethylsulfanyl)ethanethiol"),
        ("NCCCNCCS", "WR-1065"),
        ("C=CCS", "2-Propene-1-thiol"),
        ("SC(CS)(C)C", "2-methylpropane-1,2-dithiol"),
        ("SC(C(C)C)C", "3-Methyl-2-butanethiol"),
        ("C(C(S)CCO)C", "3-mercaptopentanol"),
        ("CCCC(C)SCCCO", "4-sulfanyloctan-1-ol"),  # adjusted simplified
        ("SCCCC", "butanethiol"),
        ("CC(C)(S)CCO", "3-methyl-3-sulfanylbutan-1-ol"),
        ("SC1CCCC1", "Cyclopentanethiol"),
        ("SCCCCCCCCCS", "1,9-Nonanedithiol"),
        ("SCCCCCCC", "Heptane-1-thiol"),
        ("SC(CCCCC)C", "2-Heptanethiol"),
        ("SCC(CC)C", "2-Methyl-1-butanethiol"),
        ("SC(C(C)C)CO", "xi-2-Mercapto-3-methyl-1-butanol"),
        ("OCCS", "mercaptoethanol"),
        ("SC(CCC)C", "2-Pentanethiol"),
        ("SC(CC)CS", "1,2-Butanedithiol"),
        ("SCC(CCCC)CC", "2-Ethyl-1-hexanethiol"),
        ("CC(=O)CC(C)(C)S", "4-mercapto-4-methylpentan-2-one"),
        ("CCS", "ethanethiol")
    ]
    
    for smi, name in test_cases:
        result, reason = is_alkanethiol(smi)
        print(f"SMILES: {smi:25} NAME: {name:40} --> {result}; {reason}")