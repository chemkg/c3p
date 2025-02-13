"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: Non-proteinogenic amino acid
Defined as: “Any amino acid that is not naturally encoded in the genetic code of any organism.”
This version uses several checks:
  1. It requires the presence of a carboxylic acid group.
  2. It searches for an α–carbon that is sp³ and is bound to both a carboxyl group and an [NX3] (any trivalent nitrogen),
     so that even substituted amines (N-methyl, N-hydroxy, etc.) are allowed.
  3. It excludes molecules that appear to be peptides (by looking for a peptide bond pattern).
  4. It rejects molecules whose molecular weight is far outside the typical range for single amino acids.
  5. It compares the canonical SMILES to a set of the 20 standard proteinogenic amino acids.
If the molecule has the amino acid backbone and is not one of the standard proteinogenic amino acids, 
it is classified as non–proteinogenic.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non–proteinogenic amino acid.
    
    Checks performed:
      - The molecule must parse correctly.
      - Its molecular weight should be within a typical amino acid range (roughly 75–350 Da).
      - It must have at least one carboxylic acid group.
      - It must have an alpha carbon (sp³) that is attached to both a carboxylic acid and an amino group.
        (Here, any trivalent nitrogen [NX3] is allowed so that substituted amines are accepted.)
      - Molecules containing a peptide-bond substructure are excluded.
      - Its canonical SMILES is compared to that of the 20 standard amino acids.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a non–proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Check molecular weight – typical amino acids fall roughly between 75 and 350 Da.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 75 or mw > 350:
        return False, f"Molecular weight {mw:.1f} out of range for typical amino acids"
    
    # Check for a carboxylic acid group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group detected"
    
    # Check for an alpha-carbon pattern:
    # We look for an sp3 carbon attached directly to a trivalent nitrogen ([NX3]) and to a carboxyl group.
    # This pattern should capture even N–substituted amino groups as in N,N–dihydroxydihomomethionine.
    alpha_pattern = Chem.MolFromSmarts("[C;X4]([NX3])C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(alpha_pattern):
        return False, "No typical alpha-amino acid backbone detected"
    
    # Exclude molecules that appear to be peptides.
    # Peptide bonds (an amide bond linking two alpha-carbons) are detected by this pattern.
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N[C;X4]")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Contains peptide bond; likely a peptide rather than a single amino acid"
    
    # Compute the canonical SMILES (including stereochemistry) for standardized comparison.
    input_canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    # Define canonical SMILES for the 20 standard proteinogenic amino acids.
    proteinogenic_aas = [
        "NCC(=O)O",                              # glycine
        "N[C@@H](C)C(=O)O",                       # L-alanine
        "N[C@@H](C(C)C)C(=O)O",                    # L-valine
        "N[C@@H](CC(C)C)C(=O)O",                   # L-leucine
        "N[C@@H](C[C@H](C)C)C(=O)O",               # L-isoleucine
        "O=C(O)N1CC[C@H](C1)O",                    # L-proline
        "N[C@@H](Cc1ccccc1)C(=O)O",                # L-phenylalanine
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",         # L-tryptophan
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O",              # L-tyrosine
        "N[C@@H](CO)C(=O)O",                       # L-serine
        "N[C@@H]([C@H](O)C)C(=O)O",                # L-threonine
        "N[C@@H](CS)C(=O)O",                       # L-cysteine
        "N[C@@H](CCSC)C(=O)O",                      # L-methionine
        "N[C@@H](CC(=O)O)C(=O)O",                   # L-aspartic acid
        "N[C@@H](CCC(=O)O)C(=O)O",                  # L-glutamic acid
        "N[C@@H](CC(=O)N)C(=O)O",                   # L-asparagine
        "N[C@@H](CCC(=O)N)C(=O)O",                  # L-glutamine
        "N[C@@H](CCCCN)C(=O)O",                     # L-lysine
        "N[C@@H](CCCNC(=N)N)C(=O)O",                # L-arginine
        "N[C@@H](Cc1c[nH]cn1)C(=O)O"                # L-histidine
    ]
    prot_set = set()
    for aas in proteinogenic_aas:
        aas_mol = Chem.MolFromSmiles(aas)
        if aas_mol:
            canon = Chem.MolToSmiles(aas_mol, isomericSmiles=True)
            prot_set.add(canon)
    
    if input_canonical in prot_set:
        return False, "Matches a standard proteinogenic amino acid"
    
    # If all tests pass, classify as a non-proteinogenic amino acid.
    return True, "Has amino acid functional groups and does not match any standard proteinogenic amino acid"

# For debugging or simple testing:
if __name__ == "__main__":
    # A set of test examples drawn from the outcomes.
    test_examples = [
        # True positives (non‐proteinogenic):
        ("N[C@@H](CC1=CC=C(F)C=C1)C(O)=O", "4-fluorophenyl-L-alanine"),
        ("CSCCCCCC(N(O)O)C(O)=O", "N,N-dihydroxydihomomethionine"),
        ("NC(Cc1c[nH]c2cccc(O)c12)C(O)=O", "4-hydroxytryptophan"),
        ("O(C[C@@H](N)C(O)=O)CC1=CC=CC=C1", "O-Benzyl-D-serine"),
        ("CN[C@@H](Cc1ccccc1)C(O)=O", "N-methyl-L-phenylalanine"),
        ("IC=1C=CC(C[C@H](C(O)=O)N)=CC1", "4-iodo-D-phenylalanine"),
        ("[H][C@](N)(CS[C@H](\\C=C\\C=C\\C=C/C\\C=C/CCC(O)=O)[C@@H](O)CCCC(O)=O)C(O)=O", "18-carboxy-19,20-dinor-leukotriene E4"),
        ("NC(CCCC(N)C(O)=O)C(O)CC(O)=O", "2,6-diamino-7-hydroxy-azelaic acid"),
        ("OC(=O)CNCC(O)=O", "iminodiacetic acid"),
        ("C([C@H]([C@@H](C(=O)N)O)N)(O)=O", "(3S)-3-hydroxy-L-asparagine"),
        ("NC(CC1=CC2=CC=CC=C2C=C1)C(O)=O", "3-naphthalen-2-ylalanine"),
        ("CC(NCCCC[C@H](N)C(O)=O)C(O)=O", "N(6)-(1-carboxyethyl)-L-lysine"),
        ("N[C@@H](CCCNO)C(O)=O", "N(5)-hydroxy-L-ornithine"),
        ("NCCC[C@H](N)C(O)=O", "L-ornithine"),
        ("NCCOP(O)(=O)OCC(N)C(O)=O", "serine phosphoethanolamine"),
        ("O=C(O)[C@@H](N)CCCCCCCSC", "L-pentahomomethionine"),
        ("OC([C@@H](N)CSC(=O)CCCCCCC)=O", "S-octanoyl-L-cysteine"),
        ("OC(=O)CNc1ccc(O)cc1", "N-(p-hydroxyphenyl)glycine"),
        ("CC(C)=CC\\C=C(\\C)[C@H]1CN[C@@H]([C@H]1CC(O)=O)C(O)=O", "dainic acid A"),
        ("N[C@H](Cc1c[nH]cn1)C(O)=O", "D-histidine"),
        ("N[C@H](Cc1ccc(O)c(O)c1)C(O)=O", "D-dopa"),
        ("N[C@@H](CS\\C(Cl)=C/Cl)C(O)=O", "S-(trans-1,2-dichlorovinyl)-L-cysteine"),
        ("CC(=O)CNCCCCC(N)C(O)=O", "N(6)-acetonyllysine"),
        ("N[C@@H]([C@H](O)CC(O)=O)C(O)=O", "(R)-3-hydroxy-L-glutamic acid"),
        ("NC(Cc1ccc(O)c(O)c1)C(O)=O", "dopa"),
        # Some examples that were previously (wrongly) classified:
        ("CNC(=N)NCCC[C@H](N)C([O-])=O", "N(omega)-methyl-L-argininate"),
        ("OC(=O)C(NCC(O)=O)CCCN=C(N)N", "L-Acetopine")
    ]
    for smi, name in test_examples:
        result, reason = is_non_proteinogenic_amino_acid(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nReason: {reason}\n")