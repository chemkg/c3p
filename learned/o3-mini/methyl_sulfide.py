"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl sulfide 
Definition: Any aliphatic sulfide in which at least one of the organyl groups attached 
to the sulfur is a methyl group.
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a given molecule is a methyl sulfide.
    A methyl sulfide is defined as an aliphatic sulfide (i.e. a non‐aromatic sulfur 
    atom bonded to two carbon atoms) in which at least one of the carbons attached to 
    sulfur is a methyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Make sure hydrogens are considered (for hydrogen count on carbons)
    mol = Chem.AddHs(mol)
    
    # Iterate over atoms looking for a sulfur atom that qualifies
    for atom in mol.GetAtoms():
        # Check if the atom is sulfur (atomic number 16) and non-aromatic
        if atom.GetAtomicNum() == 16 and not atom.GetIsAromatic():
            # For aliphatic sulfide we expect sulfur to have exactly 2 neighbors.
            if atom.GetDegree() != 2:
                continue
            
            neighbors = atom.GetNeighbors()
            # Ensure both attached atoms are carbons (organyl groups)
            if not all(neigh.GetAtomicNum() == 6 for neigh in neighbors):
                continue
            
            # Now check if at least one of the carbon neighbors is a methyl group.
            # We approximate a methyl group as a carbon atom that has exactly one heavy-atom neighbor (i.e. our S)
            # and three hydrogens.
            for neigh in neighbors:
                # Check that the neighbor is non-aromatic an sp3 carbon (rough approximation)
                if not neigh.GetIsAromatic() and neigh.GetAtomicNum() == 6:
                    # The number of attached hydrogens is given by GetTotalNumHs().
                    if neigh.GetTotalNumHs() == 3:
                        return True, "Found aliphatic sulfide with at least one methyl group attached."
    
    return False, "No suitable aliphatic sulfide with a methyl group attached was found."