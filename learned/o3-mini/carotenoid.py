"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python
"""
Classifies: Carotenoid
Definition:
  One of a class of tetraterpenoids (C40), formally derived from the acyclic parent,
  psi,psi-carotene by hydrogenation, dehydrogenation, cyclization, oxidation, or combination
  of these processes. This class includes carotenes, xanthophylls and certain compounds that arise
  from rearrangement of the skeleton of psi,psi-carotene or by loss of part of this structure.
  Retinoids (usually C20 derivatives) are excluded.
  
This function uses several criteria:
  1. It checks that the number of carbon atoms is sufficiently high.
  2. It searches for an extended conjugated polyene chain (using a SMARTS pattern with three adjacent C=C bonds).
  3. It optionally checks that molecular weight is in a typical range.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    A carotenoid is defined as a tetraterpenoid (derived from a C40 psi,psi-carotene skeleton)
    featuring an extended conjugated polyene chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a carotenoid, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 30:
        return False, f"Too few carbon atoms ({carbon_count} found; carotenoids generally have near 40 carbons)"
    
    # Look for an extended polyene system. Carotenoids contain long conjugated chains.
    # We require at least three consecutive C=C bonds. The SMARTS pattern "C=C-C=C-C=C"
    # seeks out a chain with three double bonds.
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No extended conjugated polyene chain found"
    
    # Optionally, check the molecular weight. Carotenoids are often >500 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a typical carotenoid"
    
    return True, f"Contains {carbon_count} carbon atoms and an extended polyene system (MW {mol_wt:.2f} Da)"

# Example usage:
if __name__ == "__main__":
    test_smiles = "CC(=O)CCCC(C)(C)C(=O)\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C"  # example carotenoid-like molecule (not necessarily valid)
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)