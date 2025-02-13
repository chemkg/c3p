"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide (gamma‐lactone with a 2‐furanone skeleton and substituted derivatives)
    
A “true” butenolide (2-furanone) in our definition should have a five‐membered lactone ring that:
  • contains exactly one ring oxygen,
  • has at least one ring carbon which bears an exocyclic C=O (i.e. a carbonyl not in the ring),
  • and that same carbon is directly connected in the ring to the ring oxygen on one side and
    on its other ring side is involved in a carbon–carbon double bond.
If a five‐membered ring in the molecule meets these conditions we call it a butenolide.
    
Note: Many substituted/fused systems complicate the picture. This algorithm attempts to minimize
false positives (molecules that incidentally match a furanone-looking unit) while still catching many of the
true butenolide variants.
    
Requires: rdkit
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide (gamma-lactone with a 2-furanone skeleton)
    based on its SMILES string.
    
    The approach:
      1. Parse the SMILES into an RDKit molecule.
      2. For each ring of size 5, check if:
           - The ring contains exactly one oxygen atom.
           - There is at least one ring carbon (atomic number 6) that has an external double bond 
             to an oxygen (carbonyl) AND in the ring it is directly connected to:
                 (a) the ring oxygen on one side, and 
                 (b) a C (atomic number 6) through a double bond.
         If such a ring is found, we call it a butenolide.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule appears to be a butenolide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each is a tuple of atom indices
    if not atom_rings:
        return False, "No rings found in molecule"

    # iterate over rings of size 5 only
    for ring in atom_rings:
        if len(ring) != 5:
            continue
        
        # Count how many ring atoms are oxygen (atomic number 8)
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue  # we require exactly one oxygen in the ring
        
        # For clarity, mark the set of indices that are in the ring.
        ring_set = set(ring)
        
        # Now search for a ring carbon that meets the following:
        #   (i) It is a carbon (atomic number 6).
        #   (ii) It has an external double bond to an oxygen (i.e. a carbonyl group).
        #   (iii) In the ring, it is directly connected to two atoms; one of these must be
        #         the ring oxygen, and the other bond (to its other ring neighbor) should be a double bond.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # looking only at carbon atoms
            
            # Find ring neighbors (neighbors whose indices are in ring_set)
            ring_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring_set]
            if len(ring_neighbors) != 2:
                continue  # in a simple ring each atom should have exactly 2 ring neighbors
            
            # Check for an exocyclic carbonyl at this carbon.
            has_carbonyl = False
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                # Skip if the neighbor is also in the ring (this is a ring bond).
                if nbr.GetIdx() in ring_set:
                    continue
                # Look for a double bond to an oxygen.
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    has_carbonyl = True
                    break
            if not has_carbonyl:
                continue  # try next ring carbon
            
            # Now check ring connectivity at this carbon.
            # We want one of its two ring neighbors to be the sole oxygen in the ring,
            # and the bond from the carbon to the other neighbor should be a C=C (double bond).
            # (Either order is acceptable.)
            nbr1 = mol.GetAtomWithIdx(ring_neighbors[0])
            nbr2 = mol.GetAtomWithIdx(ring_neighbors[1])
            # Identify if one neighbor is ring oxygen.
            if nbr1.GetAtomicNum() == 8:
                other_idx = ring_neighbors[1]
            elif nbr2.GetAtomicNum() == 8:
                other_idx = ring_neighbors[0]
            else:
                # If neither neighbor is oxygen, this carbon is not in the desired context.
                continue
            
            # Check the bond from the carbon (atom) to the "other" ring neighbor should be double.
            bond_to_other = mol.GetBondBetweenAtoms(idx, other_idx)
            if bond_to_other is None or bond_to_other.GetBondType() != Chem.BondType.DOUBLE:
                continue
            
            # If we reach here, we have found a ring carbon that:
            #   - has an external carbonyl (C=O),
            #   - is bonded in-ring to the ring oxygen on one side,
            #   - and in-ring to another carbon through a double bond.
            # Thus, we consider the ring to be a 2-furanone skeleton.
            return True, "Contains a 2-furanone skeleton typical of butenolide derivatives"
    
    # If no ring met the criteria, then we do not classify it as butenolide.
    return False, "Does not contain the required 2-furanone skeleton for butenolides"


# Example usage (remove or comment these out for production use)
if __name__ == '__main__':
    test_examples = {
        "protoanemonin": "C=C1OC(=O)C=C1",
        "Aspersclerotiorone B": "O=C1O/C(=C(/C[C@]2(OC=C(C2=O)C)C)\\C)/C(=C1)OC",
        "deslanoside": "[H][C@@]1(C[C@H](O)[C@]([H])(O[C@H]2C[C@H](O)[C@]([H])(O[C@H]3C[C@H](O)[C@]([H])(O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)[C@@H](C)O3)[C@@H](C)O2)[C@@H](C)O1)O[C@H]1CC[C@@]2(C)[C@]([H])(CC[C@]3([H])[C@]2([H])C[C@@H](O)[C@]2(C)[C@H](CC[C@]32O)C2=CC(=O)OC2",
        "5-Pentyl-3h-furan-2-one": "O1C(CCCCC)=CCC1=O",  # was false negative in previous version
        "Non-butenolide ester": "CC(=O)OC1=CC=CC=C1",  # should be false
    }
    
    for name, smi in test_examples.items():
        result, reason = is_butenolide(smi)
        print(f"{name}: {result} – {reason}")