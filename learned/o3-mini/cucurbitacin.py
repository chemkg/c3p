"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin-type compounds (tetracyclic triterpenoids derived from cucurbitane)

Heuristic notes:
1. Cucurbitacins have a tetracyclic (at least 4 rings) core.
2. They generally have a high carbon count (typically around 30 in the core skeleton, though glycosylation increases this count).
3. They are oxygenated and commonly display a conjugated enone (α,β-unsaturated carbonyl) motif.
4. The enone motif here is accepted if a carbonyl carbon (C=O) is attached by a single bond to an “alpha” carbon that is in turn doubly bonded to another carbon.
5. These criteria are heuristic and may not cover every edge‐case.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def has_enone(mol):
    """
    Checks if the molecule contains an enone motif:
    an α,β-unsaturated carbonyl. Here we search for a pattern:
       C2=C3 and C3‐C(=O) (or equivalently, a carbonyl carbon bonded via a single bond
       to a carbon that is involved in a C=C double bond).
    """
    # Loop over atoms; for each carbonyl carbon (has a C=O) search for a neighbor with a double bond.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Look for a double bond to oxygen
        has_carbonyl = False
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                has_carbonyl = True
                break
        if not has_carbonyl:
            continue  # this carbon is not a carbonyl center
        # For each carbon neighbor (alpha carbon) that is attached via a single bond...
        for bond in atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            alpha = bond.GetOtherAtom(atom)
            if alpha.GetAtomicNum() != 6:
                continue
            # Now check if alpha carbon is connected by a double bond to a beta carbon
            for bond2 in alpha.GetBonds():
                beta = bond2.GetOtherAtom(alpha)
                if beta.GetIdx() == atom.GetIdx():
                    continue
                if beta.GetAtomicNum() != 6:
                    continue
                if bond2.GetBondType() == Chem.BondType.DOUBLE:
                    return True
    return False

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin-type compound based on its SMILES string.
    
    Heuristics used:
      - Molecule must be parsed correctly.
      - It should contain at least 4 rings (the tetracyclic core; note that glycosylated derivatives may add extra rings).
      - The carbon count should be high (heuristically >=20 to allow for a 30-carbon core or its glycosylated variants).
      - The molecular weight should be above ~300 Da.
      - The compound must have an enone motif (an α,β-unsaturated carbonyl), even if the motif is integrated in a ring.
      - At least one oxygen atom is required.
      
    Returns:
        (bool, str): A tuple with True (and a reason) if the molecule meets the criteria, or False (with explanation) otherwise.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count rings.
    ring_info = mol.GetRingInfo()
    total_rings = ring_info.NumRings()
    if total_rings < 4:
        return False, f"Only {total_rings} rings detected; cucurbitacins typically are tetracyclic (>= 4 rings)"
    
    # Count carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Only {carbon_count} carbon atoms detected; too few for a typical cucurbitacin core"
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a cucurbitacin derivative"
    
    # Check for enone motif using our custom routine.
    if not has_enone(mol):
        return False, "No conjugated enone motif (α,β-unsaturated carbonyl) detected, which is common in cucurbitacins"
    
    # Count oxygen atoms.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "No oxygen atoms detected; cucurbitacins are oxygenated compounds"
    
    reason = ("Molecule has {0} rings, {1} carbons, molecular weight {2:.1f} Da, "
              "contains an enone motif, and is oxygenated, consistent with cucurbitacin-type compounds."
             ).format(total_rings, carbon_count, mol_wt)
    
    return True, reason

# Example usage (for testing):
if __name__ == '__main__':
    # Example: cucurbitacin I SMILES
    test_smiles = "CC(C)(O)\\C=C\\C(=O)[C@](C)(O)[C@H]1[C@H](O)C[C@@]2(C)[C@@H]3CC=C4[C@@H](C=C(O)C(=O)C4(C)C)[C@]3(C)C(=O)C[C@]12C"
    result, explanation = is_cucurbitacin(test_smiles)
    print(result)
    print(explanation)