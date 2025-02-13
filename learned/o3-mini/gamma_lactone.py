"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python
"""
Classifies: gamma-lactone
Definition: A gamma-lactone is a lactone having a five-membered ring.
A gamma-lactone (or furanone) has a five-membered ring containing four carbon atoms and one oxygen atom.
Within that ring one of the carbon atoms (the carbonyl carbon) should be connected to:
  - a ring oxygen via a single bond, and
  - an exocyclic oxygen via a double bond (i.e. the carbonyl group).
This program detects such a connectivity.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a cyclic ester with a five-membered ring,
    whose atoms in the ring consist of exactly four carbons and one oxygen, and
    that has one of the carboxyl carbons bearing exactly one single-bonded oxygen inside
    the ring and one double-bonded oxygen that lies outside the ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple containing True if the structure is a gamma-lactone,
                     otherwise False, along with an explanation.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all ring information as lists of atom indices
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Iterate over each ring; we are looking only at rings exactly 5 atoms long.
    for ring in ring_info:
        if len(ring) != 5:
            continue
        # Count oxygen and carbon atoms within the ring.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        ring_carbons = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        # A proper gamma-lactone ring should have exactly one oxygen and four carbons.
        if len(ring_oxygens) != 1 or len(ring_carbons) != 4:
            continue
        ring_set = set(ring)
        
        # For each carbon atom in the ring, check if it meets the lactone connectivity.
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            ring_oxygen_single = 0
            exocyclic_double = 0
            # Iterate neighbors of the candidate carbon
            for nbr in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8:
                    # Check if the oxygen is part of the ring (for the single bond)
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        if nbr.GetIdx() in ring_set:
                            ring_oxygen_single += 1
                    # Check if the oxygen is exocyclic and connected via a double bond
                    elif bond.GetBondType() == Chem.BondType.DOUBLE:
                        if nbr.GetIdx() not in ring_set:
                            exocyclic_double += 1
            # For a proper lactone the candidate carbon (the carbonyl carbon) should have
            # exactly one single-bonded ring oxygen and exactly one double-bonded exocyclic oxygen.
            if ring_oxygen_single == 1 and exocyclic_double == 1:
                return True, ("Found a five-membered ring (4C+1O) with a lactone carbonyl connectivity: "
                              "a carbon that is single-bonded to a ring oxygen and double-bonded to an exocyclic oxygen.")
    return False, "No five-membered lactone ring with the required –C(=O)O– connectivity was found"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "O1C(CCC1=O)(CCCC)CCCC",          # 5,5-Dibutyl-4,5-dihydro-2(3H)furanone (gamma-lactone)
        "O1[C@H]([C@H](CC1=O)C)CC=C(C)C",  # cis-eldanolide (gamma-lactone)
        "C=C1OC(=O)C=C1",                 # protoanemonin (gamma-lactone)
        "N[C@@H]1CCOC1=O",                # D-homoserine lactone (gamma-lactone)
        "CCCC",                          # not a lactone
    ]
    for s in test_smiles:
        result, reason = is_gamma_lactone(s)
        print(f"SMILES: {s}\n  -> Gamma-lactone? {result} | Reason: {reason}\n")