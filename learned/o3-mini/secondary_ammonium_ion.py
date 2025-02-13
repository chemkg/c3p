"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: secondary ammonium ion
Definition: An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
A protonated secondary amine (R2NH) becomes a secondary ammonium ion (R2NH2+), so we are looking for a nitrogen
with a positive formal charge, exactly two heavy (non-hydrogen) neighbors, and one hydrogen.
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string) belongs to the class of secondary ammonium ions.
    
    A secondary ammonium ion is assumed to come from protonation of a secondary amine (R2NH)
    so that, upon protonation, the nitrogen becomes R2NH2+ (two organic substituents and one hydrogen, with a +1 charge).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a secondary ammonium ion is identified, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic check to ensure the molecule is organic: it should contain at least one carbon.
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found – not an organic molecule"
    
    # Add explicit hydrogens so we can count them accurately.
    mol = Chem.AddHs(mol)
    
    # Loop over atoms to find nitrogen(s) that might be protonated secondary amines.
    for atom in mol.GetAtoms():
        # Consider only nitrogen atoms with positive formal charge.
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Get the neighbors (which now include explicit hydrogens after AddHs)
            neighbors = atom.GetNeighbors()
            heavy_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() != 1]
            hydrogen_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 1]
            
            # For a protonated secondary amine (R2NH2+),
            # we expect exactly two heavy (organic) substituents and one hydrogen attached.
            if len(heavy_neighbors) == 2 and len(hydrogen_neighbors) == 1:
                return True, "Found secondary ammonium ion: nitrogen with two organic substituents and one hydrogen (protonated secondary amine)"
    
    return False, "No protonated secondary amine (secondary ammonium ion) found"
    
# Example usage:
if __name__ == "__main__":
    test_smiles = "C[NH2+]C"  # example dimethylaminium which should classify as secondary ammonium ion.
    result, reason = is_secondary_ammonium_ion(test_smiles)
    print(result, reason)