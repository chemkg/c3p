"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide
Definition: A gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.
A butenolide is a five-membered lactone ring in which:
  - Exactly one of the five ring atoms is an oxygen (the lactone oxygen).
  - One ring carbon carries an exocyclic carbonyl group (i.e. a double bond to an oxygen that is not part of the ring).
  - At least one ring bond is unsaturated (i.e. a C=C double bond present inside the ring).
This implementation analyzes the ring information instead of relying solely on a fixed SMARTS pattern.
It should correct some of the false positives (molecules with similar substructures but not true butenolides)
and false negatives (true butenolides that do not exactly match the canonical SMARTS) from the previous version.
    
Requires: rdkit
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide (gamma-lactone with a 2-furanone skeleton)
    based on its SMILES string.
    
    Strategy:
      1. Parse the molecule.
      2. Iterate through all rings of size 5.
         - The ring must have exactly one oxygen atom (the lactone oxygen).
         - One of the ring carbons must be attached (outside the ring) to an oxygen via a double bond (the carbonyl).
         - There must be at least one double bond within the ring (indicating unsaturation).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a butenolide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()  # Tuple of tuples of atom indices for each ring
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Analyze each 5-membered ring
    for ring in ring_info:
        if len(ring) != 5:
            continue
        
        # Check for exactly one ring oxygen (lactone oxygen)
        ring_oxygen_indices = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O"]
        if len(ring_oxygen_indices) != 1:
            continue
        
        # Look for a carbon in the ring with an exocyclic double-bonded oxygen (carbonyl)
        carbonyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                # Only consider neighbors that are not in the ring
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "O":
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    # Check if bond is a double bond (we use GetBondTypeAsDouble which returns 2.0 for a double bond)
                    if bond is not None and bond.GetBondTypeAsDouble() >= 2.0:
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break

        if not carbonyl_found:
            continue

        # Check that the ring itself is unsaturated: 
        # We require at least one bond inside the ring (between consecutive ring atoms, including the closure)
        # to be a double bond. (Carbonyl double bond is exocyclic, so it won't show up as a ring bond.)
        unsaturation_found = False
        ring_atom_indices = list(ring)
        n = len(ring_atom_indices)
        # Loop over each pair of consecutive atoms (ring closure considered)
        for i in range(n):
            a_idx = ring_atom_indices[i]
            b_idx = ring_atom_indices[(i + 1) % n]
            bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
            if bond is not None and bond.GetBondTypeAsDouble() >= 2.0:
                unsaturation_found = True
                break
        
        if not unsaturation_found:
            continue
        
        # If all conditions are met for this ring, then
        return True, "Contains a 2-furanone (butenolide) five-membered lactone ring with a carbonyl and ring unsaturation"
    
    return False, "No five-membered butenolide ring found"

# Example usage (for testing purposes, remove or comment out in production)
if __name__ == '__main__':
    # Some representative examples (you can add more test cases)
    examples = {
        "protoanemonin": "C=C1OC(=O)C=C1",
        "Aspersclerotiorone B": "O=C1O/C(=C(/C[C@]2(OC=C(C2=O)C)C)\\C)/C(=C1)OC",
        "5-Pentyl-3h-furan-2-one": "O1C(CCCCC)=CCC1=O",
        "Non-butenolide ester": "CC(=O)OC1=CC=CC=C1"
    }
    
    for name, smi in examples.items():
        result, reason = is_butenolide(smi)
        print(f"{name}: {result} – {reason}")