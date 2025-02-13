"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python
"""
Classifies: gamma-lactone
Definition: A gamma-lactone is a lactone having a five-membered ring.
A gamma-lactone (or furanone) has a five-membered ring containing exactly four carbons and one oxygen.
Within that ring, one carbon (the carbonyl carbon) should be connected to:
  - a ring oxygen via a single bond, and
  - an exocyclic oxygen (i.e. not in the ring) via a double bond;
    furthermore, that exocyclic oxygen should not be further bonded.
This implementation checks for a five-membered ring with exactly one oxygen and four carbons and then
inspects each candidate carbon in the ring for the required connectivity.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a cyclic ester with a five-membered ring (4 carbons, 1 oxygen)
    where one ring carbon (the carbonyl carbon) is single-bonded to a ring oxygen and double-bonded to an exocyclic oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple containing True if the structure is classified as a gamma-lactone,
                     or False otherwise, along with an explanation.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all ring atom index lists.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # Loop over rings that are exactly 5 atoms long.
    for ring in rings:
        if len(ring) != 5:
            continue
        
        # Count the number of oxygen and carbon atoms in the ring
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        ring_carbons = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        if len(ring_oxygens) != 1 or len(ring_carbons) != 4:
            continue
        ring_set = set(ring)
        
        # For every carbon in the ring attempt to find the required lactone connectivity.
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            # Initialize counters for bonds to oxygens in two categories.
            ring_oxygen_single = 0  # Should be exactly 1: an oxygen that is a ring member via a single bond.
            exocyclic_double = 0    # Should be exactly 1: a double-bonded oxygen NOT in the ring.
            # Iterate over neighbors of this candidate carbon.
            for nbr in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # Check only for oxygen neighbors.
                if nbr.GetAtomicNum() != 8:
                    continue
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    if nbr.GetIdx() in ring_set:
                        ring_oxygen_single += 1
                elif bond.GetBondType() == Chem.BondType.DOUBLE:
                    if nbr.GetIdx() not in ring_set:
                        # Extra check: Exocyclic oxygen should be bonded only to this carbon.
                        if nbr.GetDegree() == 1:
                            exocyclic_double += 1
            # If the candidate carbon shows one ring oxygen (single bond) and one exocyclic oxygen (double bond)
            if ring_oxygen_single == 1 and exocyclic_double == 1:
                return True, ("Found a five-membered ring (4C+1O) with proper lactone carbonyl connectivity: "
                              "one ring oxygen (single bond) and one exocyclic oxygen (double bond with degree 1).")
    
    return False, "No five-membered lactone ring with the required –C(=O)O– connectivity was found"

# Example usage (only when run as a script)
if __name__ == "__main__":
    test_smiles = [
        "O1C(CCC1=O)(CCCC)CCCC",          # 5,5-Dibutyl-4,5-dihydro-2(3H)furanone --> gamma-lactone
        "O1[C@H]([C@H](CC1=O)C)CC=C(C)C",  # cis-eldanolide --> gamma-lactone
        "C=C1OC(=O)C=C1",                 # protoanemonin --> gamma-lactone
        "N[C@@H]1CCOC1=O",                # D-homoserine lactone --> gamma-lactone
        "CCCC",                          # Not a lactone
    ]
    for s in test_smiles:
        result, reason = is_gamma_lactone(s)
        print(f"SMILES: {s}\n  -> Gamma-lactone? {result} | Reason: {reason}\n")