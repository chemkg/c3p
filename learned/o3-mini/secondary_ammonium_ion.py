"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: secondary ammonium ion
Definition: An organic cation obtained by protonation of any secondary amino compound.
A protonated secondary amine (R2NH) becomes R2NH2+, meaning that the nitrogen should have 2 organic substituents (heavy atoms)
and 2 hydrogen atoms, a +1 charge, be sp3 hybridized, and ideally non‐aromatic.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string) belongs to the class of secondary ammonium ions.
    
    After protonation of a secondary amine (R2NH → R2NH2+), the nitrogen should bind to two heavy (organic) substituents
    and two hydrogen atoms. In addition, to avoid misclassifying aromatic or conjugated nitrogen atoms (e.g. imidazolium),
    we also require the nitrogen to be sp3 and non-aromatic.
    
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
    
    # Basic check: ensure the molecule contains at least one carbon (i.e. is organic)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found – not an organic molecule"
    
    # Add explicit hydrogens to count hydrogens attached to atoms properly.
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms looking for nitrogen atoms that are protonated
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Ensure the nitrogen is sp3 hybridized (as expected for protonated aliphatic amines)
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                continue
            
            # Skip aromatic nitrogen since these are typically part of oxo- or other conjugated systems.
            if atom.GetIsAromatic():
                continue
            
            # Count heavy neighbors (non-hydrogen atoms) and total attached hydrogens.
            heavy_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1)
            # GetTotalNumHs counts both implicit and explicit hydrogens attached to the atom.
            h_count = atom.GetTotalNumHs()
            
            # For a protonated secondary amine (R2NH2+), expect exactly two heavy neighbors and two hydrogens.
            if heavy_neighbors == 2 and h_count == 2:
                return True, ("Found secondary ammonium ion: protonated secondary amine with "
                              f"{heavy_neighbors} organic substituents and {h_count} hydrogens")
    
    return False, "No protonated secondary amine (secondary ammonium ion) found"
    
# Example usage:
if __name__ == "__main__":
    # Test a few cases:
    test_smiles = [
        "C[NH2+]C",  # dimethylaminium should be a secondary ammonium ion.
        "C(C[NH3+])CCC[NH2+]CC([C@@H](CO)O)=O"  # N-(D-erythrulosyl)-cadaverine: contains at least one candidate.
    ]
    for smi in test_smiles:
        result, reason = is_secondary_ammonium_ion(smi)
        print(smi, "->", result, "|", reason)