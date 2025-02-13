"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl Sulfide
Definition: Any aliphatic sulfide in which at least one of the organyl groups attached 
to the sulfur is a methyl group.
Note: In this revised version we do not exclude molecules that appear to be amino acids or peptides.
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a given molecule is a methyl sulfide.
    
    A methyl sulfide (thioether) by our definition is a molecule that contains at least one
    aliphatic (non‐aromatic) sulfur atom that is bonded to two carbon atoms and at least one 
    of these carbons is a methyl group. A methyl group is defined here as an sp3-hybridized carbon 
    that is connected to exactly one heavy (non-hydrogen) atom (the sulfur).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a methyl sulfide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are available.
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms to locate sulfur atoms.
    for atom in mol.GetAtoms():
        # Look for sulfur (atomic number 16)
        if atom.GetAtomicNum() != 16:
            continue
        # We need an aliphatic (non‐aromatic) sulfur. Skip if aromatic.
        if atom.GetIsAromatic():
            continue
        
        # We require the sulfur to be part of a thioether: it must be connected to two heavy atoms.
        # (Hydrogens are not considered heavy atoms.)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 2:
            continue
        
        # Check that both neighbors are carbons.
        if any(nbr.GetAtomicNum() != 6 for nbr in heavy_neighbors):
            continue
        
        # At least one neighbor carbon must be a methyl group.
        # We define a methyl group as a carbon (sp3-hybridized, non-aromatic) whose only heavy atom neighbor is the sulfur.
        for nbr in heavy_neighbors:
            # Skip if not sp3 or if aromatic.
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or nbr.GetIsAromatic():
                continue
            # Count heavy (non-hydrogen) neighbors of the carbon.
            nbr_heavy = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() != 1]
            # For a methyl carbon, the only heavy neighbor should be the sulfur.
            if len(nbr_heavy) == 1:
                return True, "Found an aliphatic sulfide with a methyl group attached to sulfur."
    
    return False, "No suitable aliphatic sulfide with a methyl group attached was found."

# Example usage (uncomment to test):
# test_smiles = "CCSC"  # ethyl methyl sulfide
# result, reason = is_methyl_sulfide(test_smiles)
# print(result, reason)