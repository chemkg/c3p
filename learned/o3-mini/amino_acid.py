"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: Amino acids
Definition: A carboxylic acid containing one or more amino groups.
Note:
  – The working definition here is that the molecule must have at least one carboxyl group (protonated or deprotonated)
    and at least one nitrogen atom that bears a hydrogen.
  – To avoid mis‐classifying peptides, we also count amide bonds (i.e. “N–C(=O)” connections). If the molecule has
    more than one amide bond and does not show an isolated (alpha) amino acid pattern, it is rejected.
  – This heuristic is imperfect but aims to improve on previous outcomes.
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule can be classified as an amino acid based on its SMILES string.
    Here an amino acid is any molecule that contains at least one carboxylic acid group and one or more N atoms
    with at least one hydrogen. To help reject peptides we also check for multiple amide bonds – if the molecule
    shows multiple N-C(=O) bonds and does not have an isolated alpha amino acid backbone, it is likely a peptide.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an amino acid, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for carboxylic acid group --- 
    # Two SMARTS: one for protonated –COOH and one for deprotonated –COO-
    acid_smarts1 = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_smarts2 = Chem.MolFromSmarts("[CX3](=O)[O-]")
    has_acid = mol.HasSubstructMatch(acid_smarts1) or mol.HasSubstructMatch(acid_smarts2)
    if not has_acid:
        return False, "No carboxylic acid group found"
    
    # --- Check for an amino group ---
    # Instead of excluding “N–C(=O)” (which would reject N-acyl amino acids), we simply require that
    # there is at least one nitrogen atom that has at least one hydrogen (which is true also for NH or NHR groups).
    amino_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # GetTotalNumHs counts both explicit and implicit hydrogens.
            if atom.GetTotalNumHs() > 0:
                amino_found = True
                break
    if not amino_found:
        return False, "No amino group (nitrogen with hydrogen) found"
    
    # --- Peptide check via amide bonds ---
    # We define an amide bond as the pattern [NX3]C(=O) (which will hit internal peptide bonds).
    amide_smarts = Chem.MolFromSmarts("[NX3][C](=O)")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    amide_count = len(amide_matches)
    
    # Define a SMARTS for a typical alpha amino acid backbone.
    # This covers the canonical pattern: a chiral carbon bearing an amino group and attached to a carboxyl.
    # We accept either a deprotonated or protonated carboxyl.
    alpha_smarts1 = Chem.MolFromSmarts("[C@@H](N)C(=O)[O-]")
    alpha_smarts2 = Chem.MolFromSmarts("[C@@H](N)C(=O)O")
    has_alpha_backbone = mol.HasSubstructMatch(alpha_smarts1) or mol.HasSubstructMatch(alpha_smarts2)
    
    # Heuristic: if there are multiple amide bonds (i.e. likely a peptide) and no isolated alpha backbone,
    # then reject the molecule.
    if amide_count > 1 and (not has_alpha_backbone):
        return False, "Molecule likely contains peptide bonds (multiple amide bonds) without isolated amino acid backbone"
    
    # Passed the tests; classify as an amino acid.
    return True, "Molecule contains a carboxylic acid and at least one amino group"
    
# Example usage:
if __name__ == "__main__":
    # Test with a few example SMILES strings
    tests = [
        ("OC(=O)C(N)CN1N=CC=C1", "3-(1-Pyrazolyl)-alanine"),
        ("N[C@@H](CC1=CC=C(F)C=C1)C(O)=O", "4-fluorophenyl-L-alanine"),
        ("CC(=O)Nc1c(I)cc(I)c(C(O)=O)c1I", "Acetrizoic acid"),  # N-acyl amino acid derivative
        ("O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)CCC(=O)N", "Asp-Gln-Pro (peptide)")
    ]
    
    for s, name in tests:
        result, reason = is_amino_acid(s)
        print(f"{name}: {result} -> {reason}")