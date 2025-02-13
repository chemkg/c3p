"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: Proanthocyanidin
A proanthocyanidin is defined as a flavonoid oligomer obtained by the condensation 
of two or more units of hydroxyflavans.

Heuristic criteria used here:
  1. The molecular weight must be at least 480 Da (to include dimers such as dracorubin).
  2. The molecule should contain at least 6 rings.
  3. At least two flavan-like (chroman) units are required.
     We now use a SMARTS pattern that is more general: "c1ccc2C(O)CCOc2c1"
  4. At least one pair of these flavan units must be directly connected by a bond,
     serving as a proxy for the condensation of the units.
     
Note: This approach is heuristic and may not capture every nuance of proanthocyanidin structures.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    
    A proanthocyanidin (based on our heuristic criteria) should have:
      - A molecular weight of at least 480 Da.
      - At least 6 rings.
      - At least 2 flavan-like (chroman) units as detected by a SMARTS pattern.
      - At least one bond connecting two different flavan units (indicating condensation).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as proanthocyanidin, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Check that the molecular weight is at least 480 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 480:
        return False, f"Molecular weight is {mol_wt:.1f} Da; too low for a flavonoid oligomer."
    
    # Criterion 2: The molecule should contain at least 6 rings.
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings < 6:
        return False, f"Only {num_rings} rings detected; expected at least 6 rings."
    
    # Criterion 3: Identify flavan-like units via a SMARTS pattern representing the chroman core.
    # This pattern was chosen to be more general: "c1ccc2C(O)CCOc2c1"
    flavan_smarts = "c1ccc2C(O)CCOc2c1"
    flavan_pat = Chem.MolFromSmarts(flavan_smarts)
    if flavan_pat is None:
        return False, "Error in SMARTS pattern for flavan unit."
    
    # Search for all non-overlapping matches (ignoring chirality) of the flavan substructure.
    matches = mol.GetSubstructMatches(flavan_pat, useChirality=False)
    num_flavan_units = len(matches)
    if num_flavan_units < 2:
        return False, f"Only {num_flavan_units} flavan-like unit(s) found; need at least 2 for an oligomer."
    
    # Criterion 4: Check that at least one pair of distinct flavan units is directly connected.
    connected = False
    for i in range(len(matches)):
        set_i = set(matches[i])
        for j in range(i+1, len(matches)):
            set_j = set(matches[j])
            # Skip if the two matches share any atoms.
            if set_i.intersection(set_j):
                continue
            # Check if any atom from match i is bonded to any atom from match j.
            for atom_idx_i in set_i:
                for atom_idx_j in set_j:
                    if mol.GetBondBetweenAtoms(atom_idx_i, atom_idx_j):
                        connected = True
                        break
                if connected:
                    break
            if connected:
                break
        if connected:
            break
    
    if not connected:
        return False, "Flavan-like units detected but no inter-unit bond found (lack of condensation)."
    
    reason = (
        f"Molecule has a molecular weight of {mol_wt:.1f} Da, {num_rings} rings, "
        f"and {num_flavan_units} flavan-like unit(s); at least one pair of units is connected."
    )
    return True, reason

# Example usage:
# if __name__ == '__main__':
#     sample_smiles = "O[C@H]1Cc2c(O[C@@H]1c1ccc(O)c(O)c1)cc(O)c([C@@H]1[C@@H](O)[C@H](Oc3cc(O)cc(O)c13)c1ccc(O)c(O)c1)c2O"
#     result, classification_reason = is_proanthocyanidin(sample_smiles)
#     print(result, classification_reason)