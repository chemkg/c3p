"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: Proanthocyanidin
A proanthocyanidin is defined as a flavonoid oligomer obtained by the condensation 
of two or more units of hydroxyflavans. 

Heuristic criteria used here:
  1. The molecule must have a molecular weight > 500 Da (to favor oligomers over monomers).
  2. The molecule should have multiple rings.
  3. At least two occurrences of a flavan-like unit, identified by a simplified chroman
     core SMARTS pattern "[#6]1CCOc2ccccc12", are required.
Note: This is an approximate, heuristic approach.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin, being an oligomer of hydroxyflavan (catechin-like) units, is expected to be large
    and to contain multiple flavan-like (chroman) substructures.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a proanthocyanidin, False otherwise.
        str: A reason describing the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # First criterion: check molecular weight.
    # Oligomeric proanthocyanidins (dimers or higher) are usually >500 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight is {mol_wt:.1f} Da; too low for a flavonoid oligomer."
        
    # Second criterion: check the total number of rings.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        return False, f"Only {num_rings} rings detected; expected several rings in an oligomer."
    
    # Third criterion: count the number
    # of flavan-like units using a simplified SMARTS for a chroman core.
    # A hydroxyflavan unit (like catechin) has a benzopyran (chroman) motif.
    # The SMARTS "[#6]1CCOc2ccccc12" represents a saturated 6-membered ring (with an oxygen) fused to a benzene ring.
    flavan_smarts = "[#6]1CCOc2ccccc12"
    flavan_pat = Chem.MolFromSmarts(flavan_smarts)
    if flavan_pat is None:
        return False, "Error in SMARTS pattern for flavan unit."
    
    # Use substructure search ignoring chirality to allow some variation.
    matches = mol.GetSubstructMatches(flavan_pat, useChirality=False)
    num_flavan_units = len(matches)
    if num_flavan_units < 2:
        return False, f"Only {num_flavan_units} flavan-like unit(s) found; need at least 2 for an oligomer."
        
    # If all criteria are met, classify molecule as a proanthocyanidin.
    return True, ("Molecule has a molecular weight of {:.1f} Da, {} rings, "
                  "and {} flavan-like unit match(es), consistent with a proanthocyanidin oligomer."
                 .format(mol_wt, num_rings, num_flavan_units))

# Example usage (uncomment for testing):
# if __name__ == '__main__':
#     # Example SMILES: (-)-epigallocatechin-(4beta->6)-(+)-catechin
#     sample_smiles = "O[C@H]1Cc2c(O[C@@H]1c1ccc(O)c(O)c1)cc(O)c([C@@H]1[C@@H](O)[C@H](Oc3cc(O)cc(O)c13)c2O)c2O"
#     result, reason = is_proanthocyanidin(sample_smiles)
#     print(result, reason)