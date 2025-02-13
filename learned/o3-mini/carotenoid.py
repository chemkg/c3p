"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python
"""
Classifies: Carotenoid
Definition:
  A carotenoid is one of a class of tetraterpenoids (typically near C40, though modifications may add carbons)
  featuring an extended conjugated polyene system. This class includes carotenes,
  xanthophylls and compounds arising from rearrangement or loss of portions of the psi,psi-carotene skeleton.
  Retinoids (typically C20 derivatives) are excluded.
  
Our improved classifier uses several criteria:
  1. The molecule must have at least 35 carbon atoms.
  2. Its molecular weight must be at least 500 Da.
  3. It must contain an extended conjugated polyene chain – detected either by a SMARTS 
     pattern for three consecutive C=C bonds ("C=C-C=C-C=C") or by having at least 5 non‐aromatic double bonds.
  4. The molecule should “look like” a terpenoid (mostly carbon and hydrogen). For that, we compute:
       a. Carbon fraction among all heavy atoms (atoms with atomic number > 1) must be at least 0.65.
       b. The average mass per carbon (mol_wt/carbon_count) should be in a reasonable range (≤ 18 Da).
       
If any of these tests fail, the molecule is rejected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines whether a molecule is a carotenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a carotenoid, otherwise False.
        str: A reason explaining the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Count carbon atoms
    atoms = list(mol.GetAtoms())
    carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if carbon_count < 35:
        return False, f"Too few carbon atoms ({carbon_count} found; carotenoids generally have near 40 carbons)"
        
    # Compute molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a typical carotenoid"
        
    # Check for an extended conjugated polyene chain:
    # Option 1: SMARTS pattern for three consecutive C=C bonds.
    polyene_smarts = "C=C-C=C-C=C"
    polyene_pattern = Chem.MolFromSmarts(polyene_smarts)
    has_polyene = mol.HasSubstructMatch(polyene_pattern)
    # Option 2: Count non‐aromatic double bonds.
    non_aromatic_dbl = sum(1 for bond in mol.GetBonds() 
                           if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic())
    if not (has_polyene or non_aromatic_dbl >= 5):
        return False, "No extended conjugated polyene system found"
        
    # Compute the fraction of carbon atoms among heavy atoms (non-hydrogen atoms)
    heavy_atoms = [atom for atom in atoms if atom.GetAtomicNum() > 1]
    if heavy_atoms:
        carbon_fraction = carbon_count / len(heavy_atoms)
    else:
        carbon_fraction = 0
    if carbon_fraction < 0.65:
        return False, (f"Low carbon fraction: {carbon_fraction:.2f} (expected at least 0.65 for carotenoids)")
        
    # Check average mass per carbon (mol_wt divided by carbon count).
    # Carotenoids (mostly CH compounds) typically have a value around 13–17 Da per carbon.
    avg_mass_per_carbon = mol_wt / carbon_count
    if avg_mass_per_carbon > 18:
        return False, (f"Average mass per carbon too high: {avg_mass_per_carbon:.2f} Da (suggesting extra non‐terpenoid groups)")
        
    reason = (f"Contains {carbon_count} carbon atoms, {non_aromatic_dbl} non‐aromatic double bonds, "
              f"carbon fraction {carbon_fraction:.2f}, and average mass per carbon {avg_mass_per_carbon:.2f} Da "
              f"(MW {mol_wt:.2f} Da); meets carotenoid criteria.")
    return True, reason

# For testing purposes:
if __name__ == "__main__":
    # one of the example carotenoid SMILES (luteochrome)
    test_smiles = "O1C2(C1(CCCC2(C)C))\\C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C/C=C(\\C3OC4(C(C(CCC4)(C)C)=C3)C)/C)\\C"
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)