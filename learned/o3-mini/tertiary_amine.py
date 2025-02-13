"""
Classifies: CHEBI:32876 tertiary amine
"""
#!/usr/bin/env python3
"""
Classifies: A tertiary amine

A tertiary amine is defined as “a compound formally derived from ammonia
by replacing three hydrogen atoms by hydrocarbyl groups.” In our revised
heuristic we first search for candidate nitrogen atoms that are neutral,
have no explicit hydrogen atoms, and are bonded to exactly three substituents.
Then we require that (a) every one of these three neighbors is a carbon atom,
(b) not all three are aromatic (i.e. at least one is aliphatic), and
(c) at least two substituents have “ideal” hydrocarbon character – meaning
    the carbon directly attached to the nitrogen is only attached to carbon or hydrogen.
This approach is heuristic and may still mis‐classify borderline cases.
"""

from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is of the class tertiary amine based on its SMILES string.
    
    The function looks for a candidate nitrogen atom that:
      • Is neutral and has no explicit hydrogens (H0) and exactly 3 bonds.
      • Is connected only to carbon atoms.
      • Has at least one aliphatic substituent (not all substituents are aromatic).
      • Has at least 2 substituents where the carbon attached directly to the N has
        an "ideal" hydrocarbon environment (its neighbors besides the N are only C or H).
    
    Args:
       smiles (str): SMILES string representation of the molecule.
       
    Returns:
       bool: True if a candidate tertiary amine center passes all checks, False otherwise.
       str: A reason describing the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # This SMARTS finds a neutral nitrogen with no explicit H and exactly 3 neighbors.
    # (We do not require here that each neighbor is carbon; that will be checked later.)
    pattern = Chem.MolFromSmarts("[N;!$([N+]);H0;D3]")
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No tertiary amine candidate (N with H0 and D3) found"
    
    # Helper function to judge if the substituent (the carbon directly attached to N)
    # has an "ideal" hydrocarbon neighborhood. We only inspect the immediate neighbors
    # (excluding the candidate N) and require that all are carbon or hydrogen.
    def is_ideal_substituent(carbon_atom, parentN_idx):
        for nb in carbon_atom.GetNeighbors():
            if nb.GetIdx() == parentN_idx:
                continue
            if nb.GetAtomicNum() not in (6, 1):
                return False
        return True
    
    # Now, evaluate each candidate nitrogen.
    for match in matches:
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Check that all 3 neighbors are carbon.
        neighbors = n_atom.GetNeighbors()
        if any(nb.GetAtomicNum() != 6 for nb in neighbors):
            continue  # This candidate's substituents are not all hydrocarbyl.
        
        # Count how many substituents are "ideal" and how many are aliphatic.
        ideal_count = 0
        aliphatic_count = 0
        for nb in neighbors:
            if is_ideal_substituent(nb, n_idx):
                ideal_count += 1
            # Use the atom's aromatic flag: if not aromatic, assume aliphatic.
            if not nb.GetIsAromatic():
                aliphatic_count += 1

        # Require at least one of the three groups be aliphatic.
        if aliphatic_count < 1:
            continue

        # Allow the candidate if at least 2 of the 3 substituents have ideal hydrocarbon character.
        if ideal_count < 2:
            continue

        return True, f"Found tertiary amine center at atom index {n_idx}"
    
    return False, "No valid tertiary amine center found"

# Optional testing block
if __name__ == "__main__":
    test_examples = {
        "benzphetamine": "C[C@@H](Cc1ccccc1)N(C)Cc1ccccc1",
        "lumefantrine": "CCCCN(CCCC)CC(O)c1cc(Cl)cc2\\C(=C/c3ccc(Cl)cc3)c3cc(Cl)ccc3-c12",
        "(-)-lobeline": "[H][C@]1(CCC[C@]([H])(CC(=O)c2ccccc2)N1C)C[C@H](O)c1ccccc1",
        "SYBR Green I": "C1=CC=CC=2SC(=[N+](C12)C)/C=C/3\\C=C(N(C4=CC=CC=C34)C=5C=CC=CC5)N(CCCN(C)C)CCC",
        "celesticetin": "COC(C)C(NC(=O)[C@@H]1CCCN1C)[C@H]2O[C@H](SCCOC(=O)C3=CC=CC=C3O)[C@H](O)[C@@H](O)[C@H]2O",
        "N,N-dimethylethanolamine": "CN(C)CCO",
        "triethylamine": "CCN(CC)CC",
        "(-)-selegiline": "C[C@H](Cc1ccccc1)N(C)CC#C",
        # One example previously a false positive:
        "6-chloro-4-...benzopyran-2-one": "CN(CC1=CC=CC=C1F)CC2=CC(=O)OC3=CC(=C(C=C23)Cl)O",
        # Examples that were originally missed (false negatives):
        "Tri-allate": "CC(C)N(C(C)C)C(=O)SCC(Cl)=C(Cl)Cl",
        "Pebulate": "CCCN(CCC)C(=O)SCCC",
        "EPTC": "CCCN(CCC)C(=O)SCC",
    }
    for name, smi in test_examples.items():
        res, reason = is_tertiary_amine(smi)
        print(f"{name}: {res} - {reason}")