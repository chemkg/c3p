"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: Beta-lactam
A beta-lactam is defined as a lactam (cyclic amide) in which the amide bond is contained 
within a four-membered ring (azetidinone ring). This means that the ring must include at least 
one nitrogen and one carbon attached to a double-bonded oxygen.
"""

from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam must contain a four-membered ring in which the amide bond is contained,
    meaning that the ring has at least one nitrogen and one carbon that is double-bonded to an oxygen (carbonyl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as beta-lactam, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Iterate over each ring; we are interested only in 4-membered rings.
    for ring in atom_rings:
        if len(ring) != 4:
            continue  # Skip rings that are not four-membered
        
        # Get the atoms in the candidate ring
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check if the ring has at least one nitrogen
        if not any(atom.GetAtomicNum() == 7 for atom in atoms_in_ring):
            continue
        
        # Look for a carbon atom in the ring that is part of a carbonyl group.
        # We verify that this carbon is double-bonded to an oxygen that is NOT in the ring.
        carbonyl_found = False
        for atom in atoms_in_ring:
            if atom.GetAtomicNum() == 6:  # Only consider carbon atoms for carbonyl
                for bond in atom.GetBonds():
                    # Check for a double bond to oxygen
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        neighbor = bond.GetOtherAtom(atom)
                        # Check if this neighbor is oxygen and is not in the current ring.
                        if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring:
                            carbonyl_found = True
                            break
            if carbonyl_found:
                break
        
        # If both a nitrogen and a carbonyl on a carbon in the ring are found, it's a beta-lactam.
        if carbonyl_found:
            return True, "Found beta-lactam ring: 4-membered ring with amide functionality (contains a nitrogen and a carbonyl carbon)"
    
    # No qualifying beta-lactam ring found
    return False, "No beta-lactam ring (4-membered lactam) found"