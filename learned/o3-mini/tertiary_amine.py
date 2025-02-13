"""
Classifies: CHEBI:32876 tertiary amine
"""
#!/usr/bin/env python3
"""
Classifies: A tertiary amine
A tertiary amine is defined as a compound formally derived from ammonia by replacing
three hydrogen atoms by hydrocarbyl groups.
Our approach is to first find any neutral nitrogen atoms with no directly attached hydrogens 
and exactly three neighbors that are carbon (using a SMARTS pattern). We then further require:
  • That the N is not aromatic.
  • That none of its three C neighbors is directly involved in a carbonyl (C=O) bond.
  • That not all three substituents are aromatic (i.e. at least one substituent should be aliphatic),
    consistent with the usual notion that a tertiary amine is not simply an N-(aryl)3.
If one candidate passes all tests we classify the molecule as a tertiary amine.
Note: This is a heuristic filter that may mis‐classify edge‐cases.
"""

from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is of the class tertiary amine based on its SMILES string.
    A tertiary amine is defined here as a neutral nitrogen that has no explicit hydrogen atoms,
    is connected to exactly three carbon substituents – each presumed to come from a hydrocarbyl group.
    Additional post-filtering is applied to reduce false positives. 
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule qualifies as a tertiary amine,
             False otherwise.
       str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Our SMARTS pattern finds a nitrogen atom that:
    #   - Is not positively charged (! $([N+]))
    #   - Has no attached explicit hydrogen atoms (H0)
    #   - Has exactly 3 neighbors (D3)
    #   - And each of its three substituents is a carbon atom ([$([#6])])
    tertiary_amine_smarts = "[N;!$([N+]);H0;D3]([$([#6])])([$([#6])])([$([#6])])"
    pattern = Chem.MolFromSmarts(tertiary_amine_smarts)
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No tertiary amine center found"
    
    # Iterate over each candidate nitrogen atom to perform additional filtering.
    for match in matches:
        n_idx = match[0]  # The 0th atom in the match is our nitrogen.
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # 1. The nitrogen should not be aromatic.
        if n_atom.GetIsAromatic():
            continue  # Skip candidate
        
        # 2. Among its three attached atoms (should be carbons by our SMARTS pattern),
        #    check that none is directly involved in a C=O (carbonyl) double bond.
        neighbor_idxs = [nbr.GetIdx() for nbr in n_atom.GetNeighbors()]
        has_carbonyl_neighbor = False
        aromatic_substituent_count = 0
        for nbr in n_atom.GetNeighbors():
            # We already require by SMARTS that neighbor is carbon.
            if nbr.GetAtomicNum() != 6:
                has_carbonyl_neighbor = True
                break
            # Count substituents that are aromatic.
            if nbr.GetIsAromatic():
                aromatic_substituent_count += 1
            # Check bonds from the neighbor: if any double bond to oxygen is found, mark it.
            for bond in nbr.GetBonds():
                # Get the atom at the other end of the bond.
                other = bond.GetOtherAtom(nbr)
                # If the neighbor is part of a C=O, we expect a double bond to atomic number 8.
                if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                    has_carbonyl_neighbor = True
                    break
            if has_carbonyl_neighbor:
                break
        
        if has_carbonyl_neighbor:
            continue  # This candidate N is attached to a carbonyl carbon; skip.
        
        # 3. To avoid cases where the nitrogen is attached only to aromatic carbons,
        #    at least one substituent must be aliphatic (non‐aromatic).
        if aromatic_substituent_count == 3:
            continue
        
        # If we reach here, this candidate appears acceptable.
        return True, f"Found tertiary amine center at atom index {n_idx}"
    
    # If no candidate tertiary amine passes all criteria, return False.
    return False, "No valid tertiary amine center found"

# Optional: testing code (can be removed when using as a module)
if __name__ == "__main__":
    examples = {
        "benzphetamine": "C[C@@H](Cc1ccccc1)N(C)Cc1ccccc1",
        "lumefantrine": "CCCCN(CCCC)CC(O)c1cc(Cl)cc2\\C(=C/c3ccc(Cl)cc3)c3cc(Cl)ccc3-c12",
        "(-)-lobeline": "[H][C@]1(CCC[C@]([H])(CC(=O)c2ccccc2)N1C)C[C@H](O)c1ccccc1",
        "SYBR Green I": "C1=CC=CC=2SC(=[N+](C12)C)/C=C/3\\C=C(N(C4=CC=CC=C34)C=5C=CC=CC5)N(CCCN(C)C)CCC",
        "celesticetin": "COC(C)C(NC(=O)[C@@H]1CCCN1C)[C@H]2O[C@H](SCCOC(=O)C3=CC=CC=C3O)[C@H](O)[C@@H](O)[C@H]2O",
        "N,N-dimethylethanolamine": "CN(C)CCO",
        "triethylamine": "CCN(CC)CC",
        "(-)-selegiline": "C[C@H](Cc1ccccc1)N(C)CC#C",
        # A known example that was previously a false positive:
        "6-chloro-4-...benzopyran-2-one": "CN(CC1=CC=CC=C1F)CC2=CC(=O)OC3=CC(=C(C=C23)Cl)O",
    }
    
    for name, smi in examples.items():
        res, reason = is_tertiary_amine(smi)
        print(f"{name}: {res} - {reason}")