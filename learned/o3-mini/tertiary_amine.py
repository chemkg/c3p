"""
Classifies: CHEBI:32876 tertiary amine
"""
#!/usr/bin/env python3
"""
Classifies: A tertiary amine
A tertiary amine is defined as a compound formally derived from ammonia by replacing 
three hydrogen atoms by hydrocarbyl groups. 
In our implementation we require that the molecule contains at least one nitrogen atom that 
is neutral, has no directly attached hydrogen atoms, is connected to exactly three atoms, 
and that each attachment is to a carbon atom.
"""

from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule contains a tertiary amine group based on its SMILES string.
    A tertiary amine is defined as a nitrogen (atomic number 7) that is neutral,
    has no directly attached hydrogens and is connected to three hydrocarbyl groups.
    To capture stereochemistry, rings, and other complications we use a SMARTS pattern.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule contains a tertiary amine group, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a tertiary amine.
    # This pattern looks for a nitrogen (atomic number 7) that:
    #   - Has no attached H atoms (H0)
    #   - Has exactly three neighbors (D3)
    #   - Is not charged (we do not match [N+])
    #   - Each substituent is required to be a carbon atom by using recursive SMARTS [$([#6])].
    # The pattern: [N;!$([N+]);H0;D3]([$([#6])])([$([#6])])([$([#6])])
    tertiary_amine_smarts = "[N;!$([N+]);H0;D3]([$([#6])])([$([#6])])([$([#6])])"
    pattern = Chem.MolFromSmarts(tertiary_amine_smarts)
    
    # Get all matches of the tertiary amine pattern.
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        # Return True with the index of the first matching nitrogen atom.
        tidx = matches[0][0]
        return True, f"Found tertiary amine center at atom index {tidx}"
    
    return False, "No tertiary amine center found"

# Some testing code (this section can be removed when using this as a module)
if __name__ == "__main__":
    # List of examples (both known positives and negatives) from the outcome summary.
    examples = {
        "benzphetamine": "C[C@@H](Cc1ccccc1)N(C)Cc1ccccc1",
        "lumefantrine": "CCCCN(CCCC)CC(O)c1cc(Cl)cc2\\C(=C/c3ccc(Cl)cc3)c3cc(Cl)ccc3-c12",
        "(-)-lobeline": "[H][C@]1(CCC[C@]([H])(CC(=O)c2ccccc2)N1C)C[C@H](O)c1ccccc1",
        "SYBR Green I": "C1=CC=CC=2SC(=[N+](C12)C)/C=C/3\\C=C(N(C4=CC=CC=C34)C=5C=CC=CC5)N(CCCN(C)C)CCC",
        "celesticetin": "COC(C)C(NC(=O)[C@@H]1CCCN1C)[C@H]2O[C@H](SCCOC(=O)C3=CC=CC=C3O)[C@H](O)[C@@H](O)[C@H]2O",
        "triethylamine": "CCN(CC)CC",
        "N,N-dimethylethanolamine": "CN(C)CCO",
        "(-)-selegiline": "C[C@H](Cc1ccccc1)N(C)CC#C",
        # A known negative (or false positive in previous attempt)
        "6-chloro-4-[[(2-fluorophenyl)methyl-methylamino]methyl]-7-hydroxy-1-benzopyran-2-one": "CN(CC1=CC=CC=C1F)CC2=CC(=O)OC3=CC(=C(C=C23)Cl)O",
    }
    
    for name, smi in examples.items():
        result, reason = is_tertiary_amine(smi)
        print(f"Example '{name}': {result} - {reason}")