"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: Proanthocyanidin 
A proanthocyanidin is defined as a flavonoid oligomer 
obtained by the condensation of two or more units of hydroxyflavans.
This program uses a few heuristic criteria:
  1. The molecular weight must be at least 500 Da (monomers are much smaller).
  2. The molecule should contain several ring systems.
  3. At least two substructure matches of a simplified flavan (catechin) unit are required.
Note: This is a heuristic approach and may not catch all edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is defined as an oligomer (dimer or more) 
    of hydroxyflavan units.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # First check: molecular weight.
    # Proanthocyanidins (being dimers or larger) are typically > 500 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight is {mol_wt:.1f} Da; too low for a flavonoid oligomer"
        
    # Second check: count the number of rings.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        # Most flavan units (catechin/epicatechin) contribute three rings each.
        return False, f"Only {num_rings} rings detected; expected multiple fused and separate rings for an oligomer"
    
    # Third check: count the number of flavan (catechin)–like units.
    # Here we use a simplified SMARTS pattern for a catechin (flavan-3-ol) unit.
    # This pattern looks for a bicyclic unit having an ether linkage and two hydroxyls on one aromatic ring.
    # Note: This pattern is not perfect and is intended as a heuristic.
    flavan_smarts = "O[C@@H]1Cc2c(O)cc(O)c2[C@H](O)C1"
    flavan_pat = Chem.MolFromSmarts(flavan_smarts)
    if flavan_pat is None:
        return False, "Error in SMARTS for flavan unit"
    
    # Perform substructure search with chirality turned off to allow slight mismatches.
    matches = mol.GetSubstructMatches(flavan_pat, useChirality=False)
    num_flavan_units = len(matches)
    if num_flavan_units < 2:
        return False, f"Only {num_flavan_units} flavan-like unit(s) found; need at least 2 for an oligomer"
        
    # If all checks pass, we classify the molecule as a proanthocyanidin.
    return True, ("Molecule has a molecular weight of {:.1f} Da, {} rings, "
                  "and {} flavan-like unit matches, consistent with a proanthocyanidin oligomer"
                 .format(mol_wt, num_rings, num_flavan_units))
                    
# Example usage (for testing purposes):
if __name__ == '__main__':
    # A sample SMILES from the provided list (procyanidin B6, for example)
    sample_smiles = "O[C@H]1Cc2c(O[C@@H]1c1ccc(O)c(O)c1)cc(O)c([C@H]1[C@H](O)[C@H](Oc3cc(O)cc(O)c13)c2O)c2O"
    result, reason = is_proanthocyanidin(sample_smiles)
    print(result, reason)