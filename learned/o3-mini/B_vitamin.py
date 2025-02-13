"""
Classifies: CHEBI:75769 B vitamin
"""
#!/usr/bin/env python
"""
Classifies: B vitamin compounds (vitamin B1, B2, B3, B5, B6, B7, B9, and B12)

This program employs several strategies to classify molecules as vitamin B compounds:
  1. Parse the SMILES and (if necessary) select the largest fragment.
  2. Exclude compounds that contain a typical Coenzyme A fragment.
  3. If a cobalt atom is found, immediately return True (a signature for vitamin B12 derivatives).
  4. For folate (vitamin B9) derivatives, require a pterin-like fragment PLUS either a 
     p-aminobenzoate or a glutamate fragment, and check that the molecular weight falls in the expected range.
  5. For other vitamin B compounds (B1, B2, B3, B5, B6, and B7), check that a defining SMARTS is present 
     and that the overall molecular weight is within a typical range.
  
Note:
  The SMARTS patterns and weight ranges used here are heuristic; given the chemical diversity of vitamin B compounds,
  occasional false positives or false negatives may still occur.
  
If no match is found, the function returns (False, <reason>).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is classified as a vitamin B compound based on its SMILES.
    Vitamin Bs include B1 (thiamine), B2 (riboflavin), B3 (nicotinic acid),
    B5 (pantothenic acid), B6 (pyridoxal/pyridoxamine), B7 (biotin),
    B9 (folate derivatives) and B12 (cobalamins).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if classified as a vitamin B compound, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # If multiple fragments exist (salts, counterions), select the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        # choose fragment with the most heavy atoms
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        
    # Compute molecular weight (exact)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # -----------------------------------------------------
    # Exclusion step: A common Coenzyme A fragment leads to many false positives.
    exclusion_smarts = "[SX2]CCNC(=O)"
    excl_pattern = Chem.MolFromSmarts(exclusion_smarts)
    if excl_pattern and mol.HasSubstructMatch(excl_pattern):
        return False, "Contains a Coenzyme A fragment, not a vitamin B compound"
    
    # -----------------------------------------------------
    # Vitamin B12 (cobalamins) check
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 27:
            return True, "Contains cobalt, consistent with vitamin B12 derivatives"
            
    # -----------------------------------------------------
    # Vitamin B9 (Folate) check:
    # Require both a pterin-like moiety along with either a p-aminobenzoate or a glutamate fragment.
    # and MW roughly in the range 350-700 (accounting for additional substituents).
    pterin_smarts = "n1c(=O)nc2ncnc12"  # approximate pterin fragment
    paba_smarts  = "c1ccc(cc1)C(=O)[O-]"  # p-aminobenzoate fragment (often deprotonated)
    glutamate_smarts = "[C@H](CCC(=O)[O-])C(O)=O"  # simplified glutamate fragment
    pterin_pat = Chem.MolFromSmarts(pterin_smarts)
    paba_pat  = Chem.MolFromSmarts(paba_smarts)
    glutamate_pat = Chem.MolFromSmarts(glutamate_smarts)
    
    if pterin_pat and mol.HasSubstructMatch(pterin_pat):
        if (paba_pat and mol.HasSubstructMatch(paba_pat)) or (glutamate_pat and mol.HasSubstructMatch(glutamate_pat)):
            if 350 <= mol_wt <= 700:
                return True, "Matches combined substructure patterns for folate derivatives (B9)"
    
    # -----------------------------------------------------
    # Define patterns for other vitamin B types.
    # Each tuple: (Vitamin label, SMARTS string, minimum MW, maximum MW)
    vitamin_patterns = [
        # B1 (Thiamine): The thiazolium ring is key.
        ("B1 (Thiamine)", "c1sc([n+])cn1", 250, 450),
        # B2 (Riboflavin): Isoalloxazine ring system.
        ("B2 (Riboflavin)", "c1ccc2nc3c(=O)[nH]c(=O)c3nc2c1", 350, 450),
        # B3 (Nicotinic acid): Pyridine with a carboxyl (allow either deprotonated or neutral).
        ("B3 (Nicotinic acid)", "c1ccncc1C(=O)[O-]", 100, 250),
        ("B3 (Nicotinic acid alternative)", "c1ccncc1C(=O)O", 100, 250),
        # B5 (Pantothenic acid): Aliphatic chain with terminal carboxyl.
        ("B5 (Pantothenic acid)", "CC(C)(CO)C(=O)O", 180, 260),
        # B6 (Pyridoxal): Aldehyde-functionalized pyridine.
        ("B6 (Pyridoxal)", "c1nc(CO)c(=O)c(C)cc1", 140, 250),
        # B6 (Pyridoxamine): Aminated pyridine variant.
        ("B6 (Pyridoxamine)", "c1nc(CO)cc(N)c1", 140, 250),
        # B7 (Biotin): Characteristic bicyclic ureido and thiolane substructure.
        ("B7 (Biotin)", "O=C1NC(=O)N[C@H]2CSC[C@H]12", 200, 350)
    ]
    
    # Check each vitamin pattern in turn.
    for vit_name, smarts, mw_min, mw_max in vitamin_patterns:
        try:
            pattern = Chem.MolFromSmarts(smarts)
        except Exception:
            continue
        if pattern and mol.HasSubstructMatch(pattern):
            # Allow some tolerance in weight (e.g., ±20%) to account for salt/adducts
            tol_min = mw_min * 0.8
            tol_max = mw_max * 1.2
            if tol_min <= mol_wt <= tol_max:
                return True, f"Matches substructure pattern for {vit_name}"
            else:
                # We matched a fragment but MW is off: might be a false positive.
                return False, f"Matched {vit_name} fragment but molecular weight ({mol_wt:.1f}) out of expected range ({mw_min}-{mw_max})"
                
    return False, "Does not match any known B vitamin substructure patterns"

# Example testing (uncomment the block below to run simple tests)
if __name__ == "__main__":
    test_smiles = [
        # Example vitamin B9: 6R-Tetrahydrofolic acid
        "O=C1NC(=NC=2NC[C@H](NC21)CNC3=CC=C(C=C3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)",
        # Example vitamin B9: 2-[[[4-[(2,4-diamino-6-pteridinyl)methyl-methylamino]phenyl]-oxomethyl]amino]pentanedioic acid
        "CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)NC(CCC(=O)O)C(=O)O",
        # Example vitamin B6: pyridoxal hydrochloride
        "C1(O)=C(C)N=CC(CO)=C1C([H])=O.Cl",
        # Example vitamin B1: thiamine(1+) diphosphate chloride
        "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",
        # Example vitamin B12: methylcobalamin
        "[H][C@]12[C@H](CC(N)=O)[C@@]3(C)CCC(=O)NC[C@@H](C)OP([O-])(=O)O[C@H]4[C@@H](O)[C@H](O[C@@H]4CO)n4c[n+](c5cc(C)c(C)cc45)[Co-3]456(C)N1C3=C(C)C1=[N+]4C(=CC3=[N+]5C(=C(C)C4=[N+]6[C@]2(C)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O",
    ]
    
    for smi in test_smiles:
        result, reason = is_B_vitamin(smi)
        print(f"SMILES: {smi}\nIs vitamin B? {result}\nReason: {reason}\n")