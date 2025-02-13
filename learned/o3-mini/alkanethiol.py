"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: Alkanethiol
An alkanethiol is defined as a compound in which a sulfanyl group (-SH)
is attached to an alkyl group (i.e., the sulfur is bonded to at least one sp3 carbon 
that is in a relatively simple, non-peptidic, non-ionic environment).
The approach checks for:
  • a thiol substructure defined by the SMARTS "[#16;X2&H1]-[C]"
  • that the carbon attached to S is sp3 (saturated)
  • that the molecule does not have formal charges on non-counterion atoms (to avoid ionic species)
  • that the molecule does not contain amide bonds (indicative of peptides/proteins)
  • that the molecule is not “overly complex” (e.g. many rings, which may indicate sugars, peptides, etc.)
Examples of valid alkanethiols include ethanethiol (CCS) and 1-hexanethiol (SCCCCCC).
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines whether a molecule is an alkanethiol by checking if it contains
    a thiol functional group (-SH) attached to an sp3 alkyl carbon in a relatively simple 
    (non-peptidic, non-ionic, low-complexity) environment.
    
    The function performs the following steps:
      1. Parses the SMILES string and adds explicit hydrogens.
      2. Checks that the overall molecule (ignoring obvious counterions such as Cl, Br, etc.)
         does not carry formal charges.
      3. Searches for the thiol substructure using the SMARTS "[#16;X2&H1]-[C]".
         For each found match, we verify that the carbon attached to sulfur is sp3 hybridized.
      4. Checks if there are any amide bonds (SMARTS "[NX3][CX3](=O)") which would indicate 
         a peptide or protein derivative.
      5. In addition, if the molecule is “large” (has many heavy atoms) and contains multiple rings,
         it is likely too complex to be a simple alkanethiol.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as a (simple) alkanethiol, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to ensure the –SH is explicitly represented.
    mol = Chem.AddHs(mol)
    
    # Exclude molecules carrying formal charges (except for common counter ions like halides).
    # We loop over atoms and if any atom (not a halogen) has a nonzero charge,
    # we assume the molecule is ionic/complex.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ["Cl", "Br", "F", "I"]:
            if atom.GetFormalCharge() != 0:
                return False, "Molecule has non-zero formal charges (excluding typical counterions)"
    
    # Define a SMARTS pattern for a thiol: sulfur with two bonds and one hydrogen bonded to any carbon.
    # This covers an –SH group where S is bonded to an alkyl (saturated) carbon.
    thiol_smarts = "[#16;X2&H1]-[C]"
    thiol_pattern = Chem.MolFromSmarts(thiol_smarts)
    if thiol_pattern is None:
        return False, "Invalid thiol SMARTS pattern"
    
    # Find all substructure matches for the thiol pattern.
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)
    if not thiol_matches:
        return False, "No thiol (-SH) group attached to a carbon found"
    
    # Validate that at least one of the matches has the carbon in an sp3 (saturated) environment.
    valid_thiol = False
    for match in thiol_matches:
        # match is a tuple: (S_index, C_index)
        s_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        # Check that the carbon is sp3 hybridized.
        if c_atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            valid_thiol = True
            break
    if not valid_thiol:
        return False, "Thiol group found but not attached to an sp3 (alkyl) carbon"
    
    # Exclude compounds that contain amide bonds; such bonds (pattern [NX3][CX3](=O)) are common in peptides.
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Molecule contains amide bonds, suggesting a peptide or protein derivative"
    
    # Use ring count as a proxy for complexity.
    # If the molecule has more than one ring and many heavy atoms (e.g., >20), it is likely not a simple alkanethiol.
    ring_count = Chem.GetSSSR(mol)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if ring_count > 1 and heavy_atom_count > 20:
        return False, "Molecule is too complex (multiple rings and many heavy atoms) to be a simple alkanethiol"
    
    return True, "Contains a thiol (-SH) attached to an sp3 alkyl carbon in a simple molecular context"


# Example test cases:
if __name__ == "__main__":
    test_cases = [
        ("Cl.SCCN", "Cysteamine hydrochloride"),
        ("CCS", "(ethylsulfanyl)methanethiol"),
        ("SCCCCCCS", "1,6-Hexanedithiol"),
        ("SC(C(C)C=O)CC", "3-Mercapto-2-methylpentanal"),
        ("SCCCCCC", "1-Hexanethiol"),
        ("S\\C=C(\\CC)/C", "2-Methyl-1-propenethiol"),
        ("C[C@H](O)[C@@H](S)P(O)(O)=O", "[(1R,2S)-2-hydroxy-1-thiopropyl]phosphonic acid"),
        # A few of the false positives from peptides or complex sugars:
        ("C[C@H](O)[C@@H](S)P([O-])([O-])=O", "Phosphonate derivative (charged)"),
        ("SCC(NC(=O)C(N)C(O)C)C(O)=O", "Threoninyl-Cysteine")
    ]
    for smi, name in test_cases:
        result, reason = is_alkanethiol(smi)
        print(f"SMILES: {smi}   NAME: {name}   --> {result}; {reason}")