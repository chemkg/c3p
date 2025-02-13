"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin compounds (vitamin B1, B2, B3, B5, B6, B7, B9, and B12)
Note: Vitamin B forms are structurally diverse. This implementation uses substructure
patterns designed to catch key features of some known B vitamins. These include a
check for cobalt (for B12) and a set of simplified SMARTS patterns for the others.
"""

from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a vitamin B compound based on its SMILES string.
    Vitamin B compounds include vitamin B1 (thiamine),
    B2 (riboflavin), B3 (nicotinic acid), B5 (pantothenic acid),
    B6 (pyridoxal/pyridoxamine), B7 (biotin), B9 (folates) and B12 (cobalamins).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a B vitamin, False otherwise.
        str: A reason for the classification decision.
    """

    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains a cobalt atom (atomic number: 27).
    # Vitamin B12 compounds are cobalamins and must include a cobalt ion.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 27:
            return True, "Contains cobalt, consistent with vitamin B12 derivatives"

    # List of (vitamin name, SMARTS pattern) pairs.
    # These patterns are simplified approximations.
    vitamin_patterns = [
        # Vitamin B1 (Thiamine): Contains a thiazolium ring typically with a methylated sulfur and a positively charged nitrogen.
        ("B1 (Thiamine)", "c1sc([n+])cn1"),
        # Vitamin B2 (Riboflavin): Contains the fused isoalloxazine core.
        ("B2 (Riboflavin)", "c1nc2c(c(=O)[nH]c(=O)c2[nH])n1"),
        # Vitamin B3 (Nicotinic acid): Contains a pyridine ring with a carboxyl group.
        ("B3 (Nicotinic acid)", "c1ccncc1C(=O)O"),
        # Vitamin B5 (Pantothenic acid): Contains a branched chain with a terminal carboxyl; note this is simplified.
        ("B5 (Pantothenic acid)", "CC(C)(CO)C(=O)O"),
        # Vitamin B6 (Pyridoxal / Pyridoxamine): Look for a substituted pyridine ring with an aldehyde or hydroxymethyl.
        ("B6 (Pyridoxal/Pyridoxamine)", "c1nc(CO)c(=O)c(C)cc1"),
        # Vitamin B7 (Biotin): Contains a thiophane (five‐membered ring with sulfur) fused to a ureido ring.
        ("B7 (Biotin)", "O=C1NC(=O)N[C@H]2CS[C@H]12"),
        # Vitamin B9 (Folate): Contains a pterin moiety. Here the pattern is highly simplified.
        ("B9 (Folates)", "n1cnc2c1ncnc2")
    ]

    # Try matching each vitamin substructure pattern.
    for name, smarts in vitamin_patterns:
        try:
            pattern = Chem.MolFromSmarts(smarts)
        except Exception as e:
            # If conversion fails, skip this pattern.
            continue
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, f"Matches substructure pattern for {name}"

    # If none of the known vitamin patterns (or cobalt for B12) match the molecule:
    return False, "Does not match any known B vitamin substructure patterns"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "O=C1NC(=NC=2NC[C@H](NC21)CNC3=CC=C(C=C3)C(=O)N[C@@H](CCC(O)=O)C(O)=O",  # 6R-Tetrahydrofolic acid (B9)
        "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",               # thiamine(1+) diphosphate chloride (B1)
        "[H][C@]12[C@H](CC(N)=O)[C@@]3(C)CCC(=O)NC[C@@H](C)OP([O-])(=O)O[C@H]4[...]",
        # Here you can add other test SMILES from the list (truncated for brevity)
    ]
    
    for smi in test_smiles:
        result, reason = is_B_vitamin(smi)
        print(f"SMILES: {smi}\nIs vitamin B? {result}  Reason: {reason}\n")