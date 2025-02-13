"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin compounds (vitamin B1, B2, B3, B5, B6, B7, B9, and B12)
Improved: Uses a hierarchical approach with exclusion filters and combined substructure
checks (especially for folate, B9) as well as a direct check for cobalt (for B12).
Note: Vitamin B compounds are structurally diverse. Our patterns are approximate.
"""

from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is classified as a vitamin B compound based on its SMILES.
    Vitamin Bs include B1 (thiamine), B2 (riboflavin), B3 (nicotinic acid), B5 (pantothenic acid),
    B6 (pyridoxal/pyridoxamine), B7 (biotin), B9 (folates) and B12 (cobalamins).
    
    The strategy:
      1. Parse the SMILES.
      2. Exclude molecules containing fragments known for coenzyme A (since many false positives came from CoA-containing species).
      3. Immediately return True if a cobalt atom (atomic number 27) is detected (vitamin B12 derivatives).
      4. For folates (B9), require that two key fragments (a pterin-type ring and a p-aminobenzoate group) are simultaneously present.
      5. For the others, try a list of SMARTS patterns.
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if classified as a vitamin B compound, False otherwise.
        str: A reason for the classification decision.
    """
    
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -----------------------------------------------------------
    # Exclusion check: many false positives came from Coenzyme A fragments.
    # Here we use a simplified SMARTS pattern for a CoA-like fragment.
    exclusion_smarts = "[SX2]CCNC(=O)"  # rough approximation: sulfur linked to CCNC(=O) fragments
    excl_pattern = Chem.MolFromSmarts(exclusion_smarts)
    if excl_pattern and mol.HasSubstructMatch(excl_pattern):
        return False, "Contains a coenzyme A fragment, not a vitamin B compound"
    
    # -----------------------------------------------------------
    # B12 check: Vitamin B12 (cobalamins) contain cobalt (atomic number 27).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 27:
            return True, "Contains cobalt, consistent with vitamin B12 derivatives"
    
    # -----------------------------------------------------------
    # B9 (Folate) check: Folate compounds have two distinct fragments.
    # One is the pterin or pteridine-like moiety and the other is the p-aminobenzoate.
    pterin_smarts = "[nH]c(=O)nc1ncnc1"  # simplified pattern for a pterin fragment
    paba_smarts = "c1ccc(cc1)C(=O)O"         # pattern for a para-substituted benzoate
    pterin_pattern = Chem.MolFromSmarts(pterin_smarts)
    paba_pattern = Chem.MolFromSmarts(paba_smarts)
    if pterin_pattern and paba_pattern:
        if mol.HasSubstructMatch(pterin_pattern) and mol.HasSubstructMatch(paba_pattern):
            return True, "Matches combined substructure patterns for folate (B9)"
    
    # -----------------------------------------------------------
    # Define a list of vitamin B substructure patterns (other than B9 and B12)
    vitamin_patterns = [
        # Vitamin B1 (Thiamine): The thiazolium ring is characteristic.
        ("B1 (Thiamine)", "c1sc([n+])cn1"),
        # Vitamin B2 (Riboflavin): Contains a fused isoalloxazine ring.
        ("B2 (Riboflavin)", "c1nc2c(c(=O)[nH]c(=O)c2[nH])n1"),
        # Vitamin B3 (Nicotinic acid): A pyridine ring with a carboxyl group.
        ("B3 (Nicotinic acid)", "O=C(O)c1ccncc1"),
        # Vitamin B5 (Pantothenic acid): A branched carbon chain with a terminal carboxyl.
        ("B5 (Pantothenic acid)", "CC(C)(CO)C(=O)O"),
        # Vitamin B6 (Pyridoxal/Pyridoxamine): A substituted pyridine ring with an aldehyde or alcohol.
        ("B6 (Pyridoxal/Pyridoxamine)", "c1nc(CO)c(=O)c(C)cc1"),
        # Vitamin B7 (Biotin): Features a fused ureido and thiolane (sulfur-containing) ring system.
        ("B7 (Biotin)", "O=C1NC(=O)N[C@H]2CS[C@H]12")
    ]
    
    # Iterate over the vitamin patterns
    for vit_name, smarts in vitamin_patterns:
        try:
            pattern = Chem.MolFromSmarts(smarts)
        except Exception:
            continue
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, f"Matches substructure pattern for {vit_name}"
    
    # If none of the above criteria are met, then no known vitamin B feature was found.
    return False, "Does not match any known B vitamin substructure patterns"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        # Example: 6R-Tetrahydrofolic acid (B9)
        "O=C1NC(=NC=2NC[C@H](NC21)CNC3=CC=C(C=C3)C(=O)N[C@@H](CCC(O)=O)C(O)=O",
        # Example: thiamine(1+) diphosphate chloride (B1, though might be tricky due to phosphates)
        "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",
        # Example: methylcobalamin (B12)
        "[H][C@]12[C@H](CC(N)=O)[C@@]3(C)CCC(=O)NC[C@@H](C)OP([O-])(=O)O[C@H]4[C@@H](O)[C@H](O[C@@H]4CO)n4c[n+](c5cc(C)c(C)cc45)[Co-3]456(C)N1C3=C(C)C1=[N+]4C(=CC3=[N+]5C(=C(C)C4=[N+]6[C@]2(C)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O"
    ]
    
    for smi in test_smiles:
        result, reason = is_B_vitamin(smi)
        print(f"SMILES: {smi}\nIs vitamin B? {result}\nReason: {reason}\n")