"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: Non–proteinogenic amino acid
Defined as: “Any amino acid that is not naturally encoded in the genetic code of any organism.”
This version uses several checks:
  1. The molecule must parse correctly.
  2. Its molecular weight is required to be within a range (roughly 75–700 Da).
  3. It must contain at least one carboxyl group (–C(=O)O or –C(=O)[O-]).
  4. It must contain a typical amino acid “backbone”: an sp³ (tetrahedral) carbon (the α–carbon)
     that is directly attached to (a) the carboxyl carbon (whose only other atoms are oxygen)
     and (b) an amino-type nitrogen (which is not part of an amide bond).
  5. To avoid di-/oligo-peptides, if more than one such backbone is found we reject the molecule.
  6. Finally, the canonical (isomeric) SMILES is compared to a set of the 20 standard proteinogenic amino acids.
If all tests pass, the molecule is classified as a non–proteinogenic amino acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non–proteinogenic amino acid.
    
    The following steps are taken:
      - Parse the molecule from the given SMILES.
      - Accept only molecules whose molecular weights fall roughly in the range 75–700 Da.
      - Require presence of a carboxyl group (in either acid or deprotonated form).
      - Look for a typical amino acid backbone by identifying a tetrahedral (sp³) carbon that is directly bonded
        both to:
          i) a carboxyl carbon (which in turn is bonded only to two oxygens besides the α–carbon), and
         ii) an amino nitrogen that is not itself bound to a carbonyl (i.e. not in an amide/peptide bond).
      - If more than one such backbone substructure is found, we assume the molecule is a peptide or part of a larger system.
      - Finally, if the canonical (isomeric) SMILES exactly matches one of the 20 standard proteinogenic amino acids,
        the molecule is rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a non–proteinogenic amino acid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Check molecular weight: allow from ~75 Da to ~700 Da
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 75 or mw > 700:
        return False, f"Molecular weight {mw:.1f} out of range for a typical single amino acid"
    
    # Check for the presence of a carboxyl group.
    # This SMARTS will match a carbonyl carbon with either an –OH or an –O– (i.e. protonated or deprotonated acid)
    acid_smarts = "[$([CX3](=O)[O-]),$([CX3](=O)O)]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group detected"
    
    # Look for a typical alpha-amino acid backbone.
    # We require an sp3 carbon (the α–carbon) that is bonded to:
    #   (i) an amino nitrogen that is NOT part of an amide (i.e. not followed by C=O)
    #   (ii) a carboxyl group: a carbonyl carbon that in turn is connected only to the α–carbon and two oxygens.
    # The SMARTS below uses:
    #   [#6;X4] : a tetrahedral carbon (the candidate α–carbon)
    #   ([NX3;!$(NC(=O))]) : a neighboring nitrogen that is not directly bound to a carbonyl
    #   C(=O)[O,OH] : a carboxyl group (the carboxyl carbon, with one double-bonded O and one O or OH)
    alpha_smarts = "[#6;X4]([NX3;!$(NC(=O))])C(=O)[O,OH]"
    alpha_pattern = Chem.MolFromSmarts(alpha_smarts)
    matches = mol.GetSubstructMatches(alpha_pattern)
    if len(matches) == 0:
        return False, "No typical alpha-amino acid backbone detected"
    elif len(matches) > 1:
        return False, "Contains multiple amino acid backbones; likely a peptide rather than a single amino acid"
    
    # Here we further validate the backbone match.
    # In our match tuple, we assume:
    #   match[0] = index of the candidate α–carbon,
    #   match[1] = index of the attached amine nitrogen,
    #   match[2] = index of the carboxyl carbon.
    alpha_idx, n_idx, acidC_idx = matches[0][:3]
    # Validate that the carboxyl carbon is the expected one:
    acidC = mol.GetAtomWithIdx(acidC_idx)
    # It should be bonded only to (a) the α–carbon and (b) 2 oxygens.
    oxy_neighbors = [nbr for nbr in acidC.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxy_neighbors) != 2:
        return False, "Carboxyl carbon is not bonded to exactly two oxygens"
    # Optionally, we could check that at least one of these oxygens is doubly bonded.
    # (For simplicity, we assume the SMARTS already enforces the C(=O) part.)
    
    # Exclude molecules that contain an explicit peptide bond linking α–carbons.
    # A simple pattern for a peptide bond is: C(=O)N[C]
    peptide_bond_smarts = "C(=O)N[C]"
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Contains peptide bond; likely a peptide rather than a single amino acid"
    
    # Compute the canonical (isomeric) SMILES for standardized comparison.
    input_canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    # Define the canonical SMILES for the 20 standard proteinogenic amino acids.
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
    
    return True, "Has amino acid functional groups and does not match any standard proteinogenic amino acid"


# When run as a script, we can try some example SMILES strings.
if __name__ == "__main__":
    # A few examples drawn from provided outcomes:
    test_examples = [
        # True positives (non–proteinogenic):
        ("N[C@@H](CC1=CC=C(F)C=C1)C(O)=O", "4-fluorophenyl-L-alanine"),
        ("CSCCCCCC(N(O)O)C(O)=O", "N,N-dihydroxydihomomethionine"),
        ("CSCCCCCC(NO)C(O)=O", "N-hydroxytrihomomethionine"),
        ("NC(Cc1c[nH]c2cccc(O)c12)C(O)=O", "4-hydroxytryptophan"),
        ("O(C[C@@H](N)C(O)=O)CC1=CC=CC=C1", "O-Benzyl-D-serine"),
        ("CN[C@@H](Cc1ccccc1)C(O)=O", "N-methyl-L-phenylalanine"),
        ("IC=1C=CC(C[C@H](C(O)=O)N)=CC1", "4-iodo-D-phenylalanine"),
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
        # Examples previously (wrongly) classified (we expect these to be rejected):
        ("CNC(=N)NCCC[C@H](N)C([O-])=O", "N(omega)-methyl-L-argininate"),
        ("OC(=O)C(NCC(O)=O)CCCN=C(N)N", "L-Acetopine"),
        ("N[C@@H](CCC(=O)Nc1ccccc1)C(O)=O", "N(5)-phenyl-L-glutamine")
    ]
    for smi, name in test_examples:
        result, reason = is_non_proteinogenic_amino_acid(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nReason: {reason}\n")