"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: Proanthocyanidin
A proanthocyanidin is defined as a flavonoid oligomer obtained by the condensation 
of two or more units of hydroxyflavans.

We use these heuristic criteria:
  1. Molecular weight must be at least 480 Da (to include dimers like dracorubin).
  2. The molecule should contain at least 6 rings.
  3. At least two flavan-like units are required. We search for a chroman 
     (benzopyran) core using a refined SMARTS pattern: “[#6]1CC[O]Cc2ccccc12”
  4. In addition, we require that at least one pair of these flavan units 
     is directly connected (i.e. a bond exists between an atom in one unit 
     and an atom in another). This connection is our proxy for “condensation”.
     
Note: This is still an approximate, heuristic approach.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin should have a sufficient molecular weight, a rich ring system,
    and at least two flavan (chroman) units that are connected (consistent with an oligomer).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a proanthocyanidin, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Criterion 1: molecular weight must be at least 480 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 480:
        return False, f"Molecular weight is {mol_wt:.1f} Da; too low for a flavonoid oligomer."
        
    # Criterion 2: count rings. We require at least 6 rings.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 6:
        return False, f"Only {num_rings} rings detected; expected at least 6 rings."

    # Criterion 3: identify flavan-like units via a SMARTS pattern for a chroman core.
    # This refined SMARTS pattern represents a saturated 6-membered ring containing one oxygen
    # fused to an aromatic ring.
    flavan_smarts = "[#6]1CC[O]Cc2ccccc12"
    flavan_pat = Chem.MolFromSmarts(flavan_smarts)
    if flavan_pat is None:
        return False, "Error in SMARTS pattern for flavan unit."
        
    # Use substructure search (ignoring chirality) to count the flavan units.
    # Get all non-overlapping matches.
    matches = mol.GetSubstructMatches(flavan_pat, useChirality=False)
    num_flavan_units = len(matches)
    if num_flavan_units < 2:
        return False, f"Only {num_flavan_units} flavan-like unit(s) found; need at least 2 for an oligomer."

    # Criterion 4: check connectivity between at least one pair of flavan unit matches.
    # Two units are considered connected if any atom in one match is directly bonded to any atom in the other.
    connected = False
    for i in range(len(matches)):
        set_i = set(matches[i])
        for j in range(i+1, len(matches)):
            set_j = set(matches[j])
            # If the two matches overlap, skip (we want distinct units)
            if set_i.intersection(set_j):
                continue
            # Check all bonds between atoms in set_i and set_j
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
        "Molecule has a molecular weight of {:.1f} Da, {} rings, and {} flavan-like unit match(es); "
        "at least one pair of units is connected, consistent with a proanthocyanidin oligomer."
        .format(mol_wt, num_rings, num_flavan_units)
    )
    return True, reason

# Example usage (for testing, uncomment below):
# if __name__ == '__main__':
#     sample_smiles = "O[C@H]1Cc2c(O[C@@H]1c1ccc(O)c(O)c1)cc(O)c([C@@H]1[C@@H](O)[C@H](Oc3cc(O)cc(O)c13)c2O)c2O"
#     result, classification_reason = is_proanthocyanidin(sample_smiles)
#     print(result, classification_reason)