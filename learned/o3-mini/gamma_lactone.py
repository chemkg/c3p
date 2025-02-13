"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python
"""
Classifies: gamma-lactone
Definition: A lactone having a five-membered lactone ring.
A gamma-lactone (also called a furanone) has a five-membered ring consisting of four carbon atoms and one oxygen atom.
Within that ring one carbon (the carbonyl carbon) is bonded to:
  - an exocyclic oxygen via a double bond (the carbonyl oxygen),
  - and a ring oxygen via a single bond.
This program checks that a five-membered ring has exactly one oxygen and that one of its carbon atoms shows the required connectivity.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a cyclic ester with a five-membered ring,
    whose atoms consist of four carbons and one oxygen, with one carbon (the carbonyl carbon)
    connected via:
       - a single bond to an oxygen that is part of the ring, and 
       - a double bond to an oxygen not in the ring (the carbonyl group).
    
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
    
    # Get all rings in the molecule (as lists of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Iterate over each ring
    for ring in ring_info:
        # Consider only rings of size 5 (candidate gamma-lactone rings)
        if len(ring) != 5:
            continue

        # Count the number of oxygen atoms and carbon atoms in the ring
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        ring_carbons = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]

        # A typical gamma-lactone ring should have exactly one oxygen and four carbons.
        if len(ring_oxygens) != 1 or len(ring_carbons) != 4:
            continue

        # Now try to find a carbon in the ring that meets the ester connectivity:
        # It must have a single-bonded neighbor that is the ring oxygen (inside the ring)
        # AND a double-bonded oxygen (the carbonyl) that is not part of the ring.
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Only consider carbon atoms as candidates for the carbonyl carbon.
            if atom.GetAtomicNum() != 6:
                continue

            has_ring_oxygen = False
            has_exocyclic_carbonyl = False

            for nbr in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8:
                    # Check bond type: if single and neighbor is within same ring, candidate for ring oxygen.
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        if nbr.GetIdx() in ring:
                            has_ring_oxygen = True
                    # Look for a double bond from the carbonyl carbon to oxygen.
                    elif bond.GetBondType() == Chem.BondType.DOUBLE:
                        # To avoid counting an oxygen already in the ring (should be distinct),
                        # require that the neighboring oxygen is not part of this ring.
                        if nbr.GetIdx() not in ring:
                            has_exocyclic_carbonyl = True

            if has_ring_oxygen and has_exocyclic_carbonyl:
                return True, ("Found a five-membered ring (4C+1O) with a carbonyl carbon "
                              "bearing a double-bonded exocyclic oxygen and a ring oxygen via a single bond")
    
    return False, "No five-membered lactone ring (with the required –C(=O)O– connectivity) was found"


# Example usage:
if __name__ == "__main__":
    # A few example SMILES, including known gamma-lactones and non-lactones.
    test_smiles = [
        "O1C(CCC1=O)(CCCC)CCCC",         # 5,5-Dibutyl-4,5-dihydro-2(3H)furanone (gamma-lactone)
        "O1[C@H]([C@H](CC1=O)C)CC=C(C)C", # cis-eldanolide (gamma-lactone)
        "C=C1OC(=O)C=C1",                # protoanemonin (gamma-lactone)
        "N[C@@H]1CCOC1=O",               # D-homoserine lactone (gamma-lactone)
        "CCCC",                         # not a lactone
    ]
    
    for s in test_smiles:
        result, reason = is_gamma_lactone(s)
        print(f"SMILES: {s}\n  -> Gamma-lactone? {result} | Reason: {reason}\n")