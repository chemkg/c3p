"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: Non‐proteinogenic amino acids
Defined as: “Any amino‐acid that is not naturally encoded in the genetic code of any organism.”
This improved version uses several SMARTS checks:
  1. It checks that the molecule has at least one carboxylic acid group.
  2. It checks for a free amine directly bonded to an sp3–carbon that carries an acid group.
  3. It excludes molecules that contain a peptide bond (i.e. an amide linkage between two α–carbons),
     which is typical of di‐ or oligo–peptides.
  4. It compares the canonical SMILES to a set of the 20 standard proteinogenic amino acids.
If the molecule has the amino acid functional groups and does not match any canonical amino acid, 
it is considered non–proteinogenic.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non–proteinogenic amino acid.
    
    The molecule must have at least one carboxylic acid group and one free (non–amide) amine.
    We then search for an α–amino acid substructure:
      – a free primary or secondary amine directly bonded to an sp3–carbon
        that bears a carboxylic acid group.
    To avoid classifying peptides (which may have a free N–terminus) as amino acids, 
    we check for a typical peptide–bond pattern (“C(=O)N[C]”). 
    Finally, the canonical SMILES is compared to those of the 20 canonical amino acids. 
    Only if the substructure is present and no match with a canonical amino acid is found, 
    the molecule is classified as non–proteinogenic.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a non–proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a free carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group detected"

    # Check for at least one free (non–amide) amine.
    # We require a primary or secondary amine that is not part of an amide (i.e. not bonded to a carbonyl).
    # The pattern below looks for an amine ([NX3;H2,H1]) that is connected to an sp3 carbon.
    free_amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1]-[C;X4]")
    if not mol.HasSubstructMatch(free_amine_pattern):
        return False, "No free amine group detected"
    
    # To help ensure that we are looking at a single amino acid rather than a peptide, 
    # search for a peptide–bond pattern. This pattern looks for a carbonyl linked via an amide nitrogen 
    # to an sp3–carbon (which would be the α–carbon of a second residue).
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N[C;X4]")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Contains peptide bond; likely a peptide rather than a single amino acid"
    
    # Look for an α–amino acid substructure.
    # Pattern: a free amine attached to an sp3–carbon that bears a carboxylic acid.
    # We allow the α–carbon to be non–chiral so that non–proteinogenic examples (e.g. iminodiacetic acid) are included.
    alpha_pattern = Chem.MolFromSmarts("[NX3;H2,H1]-[C;X4](C(=O)[O;H,-])")
    alt_alpha_pattern = Chem.MolFromSmarts("[NX3;H2,H1]-[C](C(=O)[O;H,-])")
    alpha_matches = mol.GetSubstructMatches(alpha_pattern)
    alpha_matches += mol.GetSubstructMatches(alt_alpha_pattern)
    if not alpha_matches:
        return False, "No typical amino acid backbone detected"
    
    # Compute the canonical SMILES (with stereochemistry) for standardized comparison.
    input_canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    # Define a set of canonical SMILES for the 20 standard proteinogenic amino acids.
    proteinogenic_aas = [
        "NCC(=O)O",                              # glycine
        "N[C@@H](C)C(=O)O",                       # L-alanine
        "N[C@@H](C(C)C)C(=O)O",                    # L-valine
        "N[C@@H](CC(C)C)C(=O)O",                   # L-leucine
        "N[C@@H](C[C@H](C)C)C(=O)O",               # L-isoleucine
        "O=C(O)N1CC[C@H](C1)O",                    # L-proline (representative cyclic form)
        "N[C@@H](Cc1ccccc1)C(=O)O",                # L-phenylalanine
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",         # L-tryptophan
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O",              # L-tyrosine
        "N[C@@H](CO)C(=O)O",                       # L-serine
        "N[C@@H]([C@H](O)C)C(=O)O",                # L-threonine (one representative)
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
    prot_smiles_set = set()
    for aas in proteinogenic_aas:
        aas_mol = Chem.MolFromSmiles(aas)
        if aas_mol:
            canon = Chem.MolToSmiles(aas_mol, isomericSmiles=True)
            prot_smiles_set.add(canon)
    
    # If the canonical SMILES exactly matches a proteinogenic amino acid, then it is not non–proteinogenic.
    if input_canonical in prot_smiles_set:
        return False, "Matches a standard proteinogenic amino acid"
    
    # If all checks pass, we classify the molecule as a non–proteinogenic amino acid.
    return True, "Has amino acid functional groups and does not match any standard proteinogenic amino acid"

# For debugging or simple testing:
if __name__ == "__main__":
    test_smiles = [
        # True positives (non-proteinogenic):
        "N[C@@H](CC1=CC=C(F)C=C1)C(O)=O",    # 4-fluorophenyl-L-alanine
        "CSCCCCC(N(O)O)C(O)=O",              # N,N-dihydroxydihomomethionine
        "OC(=O)CNCC(O)=O",                  # iminodiacetic acid (should now be detected)
        # False positive check (a dipeptide, for example):
        "O=C(NCC(=O)NC(C)C(=O)O)C(O)=O"      # a dipeptide-like structure
    ]
    for smi in test_smiles:
        result, reason = is_non_proteinogenic_amino_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")