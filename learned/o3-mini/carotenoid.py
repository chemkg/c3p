"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python
"""
Classifies: Carotenoid
Definition:
  A carotenoid is one of a class of tetraterpenoids (C40) typically derived
  from psi,psi-carotene and featuring an extended conjugated polyene system.
  This class includes carotenes, xanthophylls and compounds arising from rearrangement 
  or loss of parts of the psi,psi-carotene skeleton. Retinoids (typically C20 derivatives)
  are excluded.

This improved function uses several criteria:
  1. Checks that the total number of carbons is not too low (at least 35).
  2. Verifies that the molecular weight is in a typical carotenoid range (>=500 Da).
  3. Searches for an extended conjugated polyene system by requiring either:
       a. A match to a SMARTS pattern representing three consecutive C=C bonds, or 
       b. At least 5 non‐aromatic double bonds in the molecule.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    A carotenoid is defined as a tetraterpenoid (generally near C40, though modifications may add carbons)
    with an extended conjugated polyene system.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a carotenoid, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 35:
        return False, f"Too few carbon atoms ({carbon_count} found; carotenoids generally have near 40 carbons)"
        
    # Optionally, check for a typical molecular weight (most carotenoids are >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for typical carotenoids"

    # Criterion 1: Look for an extended conjugated polyene chain.
    # We use a SMARTS pattern for three consecutive C=C bonds.
    polyene_smarts = "C=C-C=C-C=C"
    polyene_pattern = Chem.MolFromSmarts(polyene_smarts)
    has_polyene_match = mol.HasSubstructMatch(polyene_pattern)
    
    # Criterion 2: Alternatively, count the number of non‐aromatic double bonds.
    # (Some molecules might have interruptions that break the strictly contiguous chain.)
    non_aromatic_dbl = sum(1 for bond in mol.GetBonds() 
                           if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic())
    
    # We consider that an extended conjugated system is present if either:
    # - The polyene SMARTS is matched, OR
    # - There are at least 5 non‐aromatic double bonds in the molecule.
    if not (has_polyene_match or non_aromatic_dbl >= 5):
        return False, "No extended conjugated polyene system found"
    
    # If all criteria are met, classify as carotenoid.
    reason = (f"Contains {carbon_count} carbon atoms, {non_aromatic_dbl} non-aromatic double bonds, "
              f"and an extended polyene system (MW {mol_wt:.2f} Da)")
    return True, reason

# Example usage:
if __name__ == "__main__":
    # Example SMILES (one of the structures that belongs to the carotenoid class)
    test_smiles = "O1C2(C1(CCCC2(C)C)C)\\C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C/C=C(\\C3OC4(C(C(CCC4)(C)C)=C3)C)/C)\\C"
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)