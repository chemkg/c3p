"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin compounds (vitamin B1, B2, B3, B5, B6, B7, B9, and B12)
Improved version:
- Excludes Coenzyme A-like fragments.
- Immediately returns True if a cobalt atom is present (a signature of B12 derivatives).
- Applies a combined substructure check for folate (B9) derivatives using a pterin-like fragment
  and either a p-aminobenzoate or glutamate-like fragment.
- Uses several specific SMARTS patterns for the remaining B vitamins.
Note: Due to the diversity of vitamin B molecules, these patterns are still approximate.
"""

from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is classified as a vitamin B compound based on its SMILES.
    Vitamin Bs include B1 (thiamine), B2 (riboflavin), B3 (nicotinic acid),
    B5 (pantothenic acid), B6 (pyridoxal/pyridoxamine), B7 (biotin),
    B9 (folate derivatives) and B12 (cobalamins).

    The strategy:
      1. Parse the SMILES string.
      2. Exclude compounds with fragments suggestive of Coenzyme A.
      3. Immediately return True if a cobalt atom (atomic number 27) is detected.
      4. For folates (B9), require the presence of a pterin-like fragment plus either a 
         p-aminobenzoate or glutamate fragment.
      5. Check for several specific vitamin B SMARTS patterns (B1, B2, B3, B5, B6, and B7).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if classified as a vitamin B compound, False otherwise.
        str: Reason for classification.
    """
    
    # Convert SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # -----------------------------------------------------------
    # Exclusion: Remove obvious Coenzyme A fragments.
    # Many false positives were coming from molecules containing a fragment like:
    # a sulfur attached to CCNC(=O) groups.
    exclusion_smarts = "[SX2]CCNC(=O)"
    excl_pattern = Chem.MolFromSmarts(exclusion_smarts)
    if excl_pattern and mol.HasSubstructMatch(excl_pattern):
        return False, "Contains a Coenzyme A fragment, not a vitamin B compound"
        
    # -----------------------------------------------------------
    # B12 Check: Vitamin B12 derivatives (cobalamins) contain cobalt (atomic number 27).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 27:
            return True, "Contains cobalt, consistent with vitamin B12 derivatives"
    
    # -----------------------------------------------------------
    # B9 (Folate) Check: Folate compounds typically contain:
    #   - A pterin-type moiety (our pattern here is a simplified pterin-like fragment)
    #   - A p-aminobenzoate moiety OR a glutamate-like fragment commonly attached to folates.
    pterin_smarts = "[nH]c(=O)nc1ncnc1"  # approximate pterin fragment
    paba_smarts  = "c1ccc(cc1)C(=O)[O-]"   # p-aminobenzoate substructure, deprotonated carboxyl
    glutamate_smarts = "[C@H](CCC(=O)O)C(O)=O"  # simplified glutamate fragment
    pterin_pat = Chem.MolFromSmarts(pterin_smarts)
    paba_pat  = Chem.MolFromSmarts(paba_smarts)
    glutamate_pat = Chem.MolFromSmarts(glutamate_smarts)
    if pterin_pat:
        if (paba_pat and mol.HasSubstructMatch(paba_pat)) or (glutamate_pat and mol.HasSubstructMatch(glutamate_pat)):
            if mol.HasSubstructMatch(pterin_pat):
                return True, "Matches combined substructure patterns for folate derivatives (B9)"
    
    # -----------------------------------------------------------
    # Define vitamin B substructure patterns (for B1, B2, B3, B5, B6, and B7).
    # Note: These SMARTS are approximate and may have variants.
    vitamin_patterns = [
        # Vitamin B1 (Thiamine): Characteristic thiazolium ring.
        ("B1 (Thiamine)", "c1sc([n+])cn1"),
        # Vitamin B2 (Riboflavin): Fused isoalloxazine ring system.
        ("B2 (Riboflavin)", "c1ccc2nc3c(=O)[nH]c(=O)c3nc2c1"),
        # Vitamin B3 (Nicotinic acid): Pyridine ring with a carboxyl group.
        ("B3 (Nicotinic acid)", "c1ccncc1C(=O)[O-]"),
        # Vitamin B5 (Pantothenic acid): Contains a branched aliphatic chain terminating in a carboxyl.
        ("B5 (Pantothenic acid)", "CC(C)(CO)C(=O)O"),
        # Vitamin B6 (Pyridoxal): Aldehyde-functionalized pyridine derivatives.
        ("B6 (Pyridoxal)", "c1nc(CO)c(=O)c(C)cc1"),
        # Vitamin B6 (Pyridoxamine): Amino-substituted pyridine variant.
        ("B6 (Pyridoxamine)", "c1nc(CO)cc(N)c1"),
        # Vitamin B7 (Biotin): Fused ureido and thiolane ring system.
        ("B7 (Biotin)", "O=C1NC(=O)N[C@H]2CSC[C@H]12")
    ]
    
    # Iterate through the defined vitamin patterns.
    for vit_name, smarts in vitamin_patterns:
        try:
            pattern = Chem.MolFromSmarts(smarts)
        except Exception:
            continue
        if pattern and mol.HasSubstructMatch(pattern):
            return True, f"Matches substructure pattern for {vit_name}"
    
    return False, "Does not match any known B vitamin substructure patterns"

# Uncomment below for simple testing.
if __name__ == "__main__":
    test_smiles = [
        # Example: 6R-Tetrahydrofolic acid (B9)
        "O=C1NC(=NC=2NC[C@H](NC21)CNC3=CC=C(C=C3)C(=O)N[C@@H](CCC(O)=O)C(O)=O",
        # Example: thiamine(1+) diphosphate chloride (B1 example)
        "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",
        # Example: methylcobalamin (B12)
        "[H][C@]12[C@H](CC(N)=O)[C@@]3(C)CCC(=O)NC[C@@H](C)OP([O-])(=O)O[C@H]4[C@@H](O)[C@H](O[C@@H]4CO)n4c[n+](c5cc(C)c(C)cc45)[Co-3]456(C)N1C3=C(C)C1=[N+]4C(=CC3=[N+]5C(=C(C)C4=[N+]6[C@]2(C)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O"
    ]
    
    for smi in test_smiles:
        result, reason = is_B_vitamin(smi)
        print(f"SMILES: {smi}\nIs vitamin B? {result}\nReason: {reason}\n")