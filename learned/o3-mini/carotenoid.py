"""
Classifies: CHEBI:23044 carotenoid
"""
#!/usr/bin/env python
"""
Classifies: Carotenoid

Definition:
  A carotenoid is a tetraterpenoid (typically near C40, though modifications may add carbons)
  featuring an extended conjugated polyene system. This class includes carotenes and xanthophylls
  as well as compounds arising from rearrangement or partial loss of the psi,psi-carotene skeleton.
  Retinoids (typically C20) are excluded.

Improved criteria used in this classifier:
  1. The molecule must have at least 35 carbon atoms.
  2. Its molecular weight must be at least 500 Da.
  3. It must have an extended conjugated polyene system:
       - Detected by a SMARTS pattern for three consecutive C=C bonds (i.e. a chain of 3 doubles with no interruption)
       - And with at least 5 non‐aromatic double bonds overall.
  4. The molecule should “look like” a terpenoid:
       a. The fraction of carbon atoms among heavy atoms (atomic number>1) is at least 0.65.
       b. The average mass per carbon (mol_wt/carbon_count) should be ≤ 18 Da.
       
If any test fails the molecule is rejected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines whether a given molecule is a carotenoid based on its SMILES string.
    
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
    
    # Get list of atoms for further computations
    atoms = list(mol.GetAtoms())
    
    # Count carbon atoms
    carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if carbon_count < 35:
        return False, f"Too few carbon atoms ({carbon_count} found; carotenoids generally have near 40 carbons)"
    
    # Calculate the molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a typical carotenoid"
    
    # Check for an extended conjugated polyene chain.
    # We require a contiguous chain of three C=C bonds (i.e. 6 sp2 carbons connected in a row via double bonds and single bonds).
    # The SMARTS below finds a pattern: carbon double-bonded to carbon, single-bond (connecting two sp2 carbons), then repeating twice.
    polyene_smarts = "[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"
    polyene_pattern = Chem.MolFromSmarts(polyene_smarts)
    has_polyene = mol.HasSubstructMatch(polyene_pattern)
    
    # Next, count all non‐aromatic double bonds (which may be part of conjugated systems)
    non_aromatic_dbonds = sum(1 for bond in mol.GetBonds() 
                              if (bond.GetBondType() == Chem.BondType.DOUBLE) and not bond.GetIsAromatic())
    
    if not has_polyene or non_aromatic_dbonds < 5:
        return False, ("No extended conjugated polyene system detected; "
                       "need at least one contiguous chain of three C=C bonds and at least 5 non‐aromatic double bonds overall")
    
    # Compute carbon fraction among heavy atoms (atomic number > 1; i.e. non-H atoms)
    heavy_atoms = [atom for atom in atoms if atom.GetAtomicNum() > 1]
    if heavy_atoms:
        carbon_fraction = carbon_count / len(heavy_atoms)
    else:
        carbon_fraction = 0
    if carbon_fraction < 0.65:
        return False, f"Low carbon fraction: {carbon_fraction:.2f} (expected at least 0.65 for carotenoids)"
    
    # Check average mass per carbon. For mainly CH compounds typical for carotenoids, value is usually 13–18 Da.
    avg_mass_per_carbon = mol_wt / carbon_count
    if avg_mass_per_carbon > 18:
        return False, f"Average mass per carbon too high: {avg_mass_per_carbon:.2f} Da (suggesting extra non‐terpenoid groups)"
    
    reason = (f"Contains {carbon_count} carbon atoms, {non_aromatic_dbonds} non‐aromatic double bonds, "
              f"carbon fraction {carbon_fraction:.2f}, and average mass per carbon {avg_mass_per_carbon:.2f} Da "
              f"(MW {mol_wt:.2f} Da); meets carotenoid criteria.")
    return True, reason


# For testing purposes:
if __name__ == "__main__":
    # Example: (5R,5'R,6S,8'R)-Luteochrome SMILES provided in the test set
    test_smiles = "O1C2(C1(CCCC2(C)C))\\C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C/C=C(\\C3OC4(C(C(CCC4)(C)C)=C3)C)/C)\\C"
    result, explanation = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", explanation)