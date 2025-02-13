"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam Antibiotic 
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
A beta-lactam ring is a four-membered cyclic amide (azetidinone) with one nitrogen and one carbonyl carbon.
"""

from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    The test is done by analyzing all four-membered rings in the molecule. For any 4-membered ring,
    if it contains exactly one nitrogen atom and at least one carbon that has a double-bonded oxygen (i.e.
    a carbonyl group), the molecule is considered to have a beta-lactam ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a beta-lactam antibiotic, False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule. This returns a tuple of tuples with atom indices for each ring.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Iterate over all rings in the molecule
    for ring in atom_rings:
        if len(ring) == 4:  # Focus on 4-membered rings only
            n_count = 0      # Counter for nitrogen atoms in the ring
            carbonyl_found = False  # Flag if a carbon with a double-bonded oxygen is present
            
            # Check each atom within the ring
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Count the number of nitrogens
                if atom.GetAtomicNum() == 7:
                    n_count += 1
                # For carbons, check if one of their bonds is a double bond to an oxygen (carbonyl)
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # Check for double bond and that the other atom is oxygen
                        if bond.GetBondTypeAsDouble() == 2.0:
                            other_atom = bond.GetOtherAtom(atom)
                            if other_atom.GetAtomicNum() == 8:
                                # Additional check: ensure the oxygen is not also part of the ring,
                                # but even if it is, the presence of a carbonyl is significant.
                                carbonyl_found = True
                                break  # No need to check other bonds for this atom
                # Early exit if conditions are already met for this ring
                if n_count > 1 and carbonyl_found:
                    # While beta-lactam ring ideally has only one N, if there are more, this ring might not be valid.
                    # We continue checking other rings.
                    break
            
            # If exactly one nitrogen and a carbonyl were found in this 4-membered ring, classify as beta-lactam.
            if n_count == 1 and carbonyl_found:
                return True, ("Molecule contains a four-membered ring with one nitrogen and a carbonyl carbon, "
                              "which is consistent with a beta-lactam ring.")
    
    # If no 4-membered ring meeting the criteria is found, reject the molecule.
    return False, "Molecule does not contain a beta-lactam ring."

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "O=C1CN(C1)C2=CC=CC=C2"  # A simple beta-lactam derivative example
    result, reason = is_beta_lactam_antibiotic(test_smiles)
    print(result, reason)