"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary Amine
Definition: A compound formally derived from ammonia by replacing one hydrogen by a hydrocarbyl group.
For this classifier, the molecule is considered a (simple) primary amine if it satisfies ALL of the following:
  • It does NOT contain any amide bonds (a substructure matching "C(=O)N") – these usually indicate peptides/amides.
  • It does NOT contain a carboxyl group (pattern "C(=O)[OH]") – many amino acids are not considered here.
  • It does NOT contain any phosphorus atoms (which are typical of e.g. phospholipids).
  • It contains at least one nitrogen atom (atomic number 7) that (after adding explicit Hs) has exactly two hydrogen atoms 
    attached AND exactly one heavy (i.e. non‐hydrogen) neighbor.
This nitrogen is interpreted as being derived from ammonia (NH3) by a single substitution with a hydrocarbyl group – that is, a primary amine R–NH2.
Note: This approach is a heuristic (not perfect) and may not correctly classify all borderline molecules.
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a (simple) primary amine based on its SMILES string.
    The classifier follows these rules:
      - Reject molecules that contain amide bonds (SMARTS "C(=O)N"), carboxylic acid groups (SMARTS "C(=O)[OH]"), or phosphorus atoms.
      - Then, look for nitrogen atoms (after adding explicit hydrogens) that have exactly 2 hydrogen atoms attached
        and exactly 1 connected heavy (non-hydrogen) atom.
      - If at least one such nitrogen is found, the molecule is considered a primary amine.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a primary amine, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we can count H atoms properly.
    mol = Chem.AddHs(mol)
    
    # Reject if the molecule contains any amide bond.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if amide_smarts and mol.HasSubstructMatch(amide_smarts):
        return False, "Molecule contains amide bond(s) typical of peptides/amides"
    
    # Reject if the molecule contains carboxylic acid group (commonly found in amino acids)
    acid_smarts = Chem.MolFromSmarts("C(=O)[OH]")
    if acid_smarts and mol.HasSubstructMatch(acid_smarts):
        return False, "Molecule contains a carboxyl group, not a simple primary amine"
    
    # Reject if the molecule contains phosphorus (e.g. phospholipids)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus atoms, likely not a primary amine"
    
    # Now, inspect each nitrogen atom to check if it is a valid primary amine:
    # (i) atomicNum == 7,
    # (ii) exactly two hydrogens attached,
    # (iii) exactly one non-hydrogen neighbor.
    valid_primary_amine_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # skip non-nitrogen
        
        # Count total number of hydrogen atoms attached.
        hydrogen_count = atom.GetTotalNumHs()
        
        # Count number of heavy (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        heavy_count = len(heavy_neighbors)
        
        # For a valid primary amine:
        # it should be derived from NH3 by replacing ONE H by a hydrocarbyl group,
        # so it must have exactly two hydrogens (NH2) and exactly one heavy neighbor.
        if hydrogen_count == 2 and heavy_count == 1:
            valid_primary_amine_found = True
            break

    if valid_primary_amine_found:
        return True, "Molecule contains a valid primary amine (R–NH2) group (nitrogen with 2 H's and 1 heavy substituent)"
    else:
        return False, "No valid primary amine group (with exactly 2 hydrogens and 1 non-hydrogen substituent) found"

# For quick testing (remove or comment-out in production)
if __name__ == "__main__":
    test_smiles_list = [
        # True positives
        "Cc1cc(cc(c1N)S(O)(=O)=O)C(=C1\\C=CC(=N)C(=C1)S(O)(=O)=O)\\c1ccc(N)c(c1)S(O)(=O)=O",  # acid fuchsin (free acid form)
        "CC(C)(C)NC[C@H](O)c1cc(Cl)c(N)c(Cl)c1",  # (R)-clenbuterol
        "NC[C@H](CC)C",  # 2-Methylbutylamine
        "CC(C)(C)NCC(O)c1cc(Cl)c(N)c(Cl)c1",  # clenbuterol
        "NNCCc1ccccc1",  # Phenelzine
        "[H]C(=O)CCCNCCCN",  # N-(3-aminopropyl)-4-aminobutanal
        "NCCCCC(CN)CCCN",  # 4-(aminomethyl)octane-1,8-diamine
        "Nc1ccc(cc1)N=Nc1ccc(cc1)[N+]([O-])=O",  # 4-(4-nitrophenylazo)aniline
        "Nc1ccc(cc1S(O)(=O)=O)\\N=N\\c1ccc(cc1)S(O)(=O)=O",  # 4-aminoazobenzene-3,4'-disulfonic acid
        "CCCCCCCCCCCCCCCN",  # 1-aminopentadecane
        "CCC(C)CC(C)N",  # methylhexaneamine
        "Nc1cnnc(O)c1Cl",  # chloridazone-desphenyl
        "C1=CC(=C(C(=C1)O)O)CCN",  # 3-(2-aminoethyl)benzene-1,2-diol
        "Cc1cc(ccc1N=Nc1ccc2c(cc(c(N)c2c1O)S(O)(=O)=O)S(O)(=O)=O)-c1ccc(N=Nc2ccc3c(cc(c(N)c3c2O)S(O)(=O)=O)S(O)(=O)=O)c(C)c1",  # Evans blue free acid
        "C[C@@H](N)Cc1ccccc1",  # (R)-amphetamine
        "Nc1ccc2ccccc2c1",  # 2-naphthylamine
        "CN",  # methylamine
        "NC(CC(C)(C)C)(C)C",  # 2,4,4-trimethyl-2-Pentanamine
        "C[C@H](N)c1ccccc1",  # (1S)-1-phenylethanamine
        "Nc1ccc2cc3ccccc3cc2c1",  # anthracen-2-amine
        "CCCCCN",  # pentan-1-amine
        "[C@@H]1(C2=CC=CC=C2)[C@@H](N)C1",  # (1S,2R)-tranylcypromine
        "CCCCN",  # butan-1-amine
        "NCCCCCCCCCCC[C@H](CC)C",  # Medelamine B
        "NC1=C(C=CC=C1)C1=NC2=C(N1)C=CC=C2",  # 2-(1H-benzimidazol-2-yl)aniline
        # False positives (the following are expected to be rejected)
        "O(CCOCCN)CCN",  # ethylene glycol bis(2-aminoethyl) ether
        "CSCCCCCCC(N)C(O)=O",  # tetrahomomethionine
        "N[C@@H](CC(=O)c1ccccc1N)C(O)=O",  # L-kynurenine
        "C[C@H](O)[C@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O",  # L-erythro-7,8-dihydrobiopterin
        # False negatives (examples where the primary amine might be “hidden” in a complex structure)
        "CC1=CC=CC=C1CCN2CCCC(C2)CN(C)CC3=CC(=CC=C3)OC",  # N-[(3-methoxyphenyl)methyl]-N-methyl-1-[1-[2-(2-methylphenyl)ethyl]-3-piperidinyl]methanamine
    ]
    
    for smi in test_smiles_list:
        result, reason = is_primary_amine(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")