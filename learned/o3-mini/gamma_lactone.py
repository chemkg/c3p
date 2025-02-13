"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python
"""
Classifies: gamma-lactone
Definition: A lactone having a five-membered lactone ring.
A gamma-lactone contains a five-membered ring in which an ester functional group (–C(=O)O–) is embedded:
the carbonyl carbon (in the ring) is double-bonded to an oxygen (typically exocyclic)
and single-bonded to an oxygen that is part of the ring.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a lactone (cyclic ester) that contains a five-membered ring.
    
    The algorithm:
      1. Parse the SMILES string into a molecule.
      2. Iterate over all rings and focus on rings with exactly 5 atoms.
      3. For each 5-membered ring, check if it contains a carbon atom
         that is bonded to:
           - an oxygen atom (with a single bond) that is part of the ring (the ring-embedded ester oxygen)
           - and an oxygen atom with a double bond (the carbonyl oxygen, which need not be in the ring).
         If this is found in any five-membered ring the molecule is classified as a gamma-lactone.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple containing True if the structure is a gamma-lactone, otherwise False,
                     plus a reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the ring information for the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"

    # Iterate over all rings
    for ring in ring_info:
        if len(ring) == 5:  # consider only 5-membered rings
            # For each atom index in the ring
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                # Look for a carbon atom (atomic number 6)
                if atom.GetAtomicNum() != 6:
                    continue

                # Initialize flags for the required connectivity.
                has_ring_oxygen = False  # oxygen in the ring via a single bond
                has_carbonyl = False     # double bonded oxygen (carbonyl) neighbor
                # Iterate over neighbors of the candidate carbon atom
                for nbr in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    # Check for oxygen neighbor
                    if nbr.GetAtomicNum() == 8:
                        # Check bond type
                        if bond.GetBondType() == Chem.BondType.SINGLE:
                            # For lactone, the oxygen connected by single bond must be in the same ring
                            if nbr.GetIdx() in ring:
                                has_ring_oxygen = True
                        elif bond.GetBondType() == Chem.BondType.DOUBLE:
                            # Found a carbonyl oxygen (does not need to be in ring)
                            has_carbonyl = True
                # If both conditions are satisfied then we have a lactone ester motif
                if has_ring_oxygen and has_carbonyl:
                    return True, "Found a five-membered ring containing the ester (lactone) motif (–C(=O)O–)"
    
    return False, "No five-membered lactone ring (with the required –C(=O)O– functionality) was found"
    
# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "O1C(CCC1=O)(CCCC)CCCC",  # 5,5-Dibutyl-4,5-dihydro-2(3H)furanone
        "C=C1OC(=O)C=C1",         # protoanemonin
        "N[C@@H]1CCOC1=O",        # D-homoserine lactone
        "CCCC",                   # not a lactone
    ]
    
    for s in test_smiles:
        result, reason = is_gamma_lactone(s)
        print(f"SMILES: {s}\n  -> Gamma-lactone? {result} | Reason: {reason}")