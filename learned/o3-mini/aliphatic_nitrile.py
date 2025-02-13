"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: Aliphatic Nitrile
Definition: Any nitrile derived from an aliphatic compound.
This program checks for the presence of a nitrile (C≡N) group and then
verifies that the nitrile group is attached to an aliphatic (non‐aromatic)
environment. In our improved approach we require that the nitrile carbon is non‐aromatic,
and that the single substituent (the atom bound to the nitrile carbon other than the nitrogen)
and its immediate neighbors do not show aromaticity.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    We require that at least one nitrile group (C≡N) is present whose nitrile carbon
    is not aromatic and is attached to a substituent whose immediate environment is free
    from aromatic atoms. This helps to distinguish nitriles derived from aliphatic,
    non‐aromatic precursors from those attached directly to aromatic systems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an aliphatic nitrile, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We'll use a SMARTS pattern to find nitrile groups.
    # Our pattern "[C]#N" will match any carbon triple bonded to nitrogen.
    # Later we will explicitly check that the nitrile carbon is non-aromatic.
    nitrile_pattern = Chem.MolFromSmarts("[C]#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group found in the molecule"
    
    ringinfo = mol.GetRingInfo()  # to check ring aromaticity

    # Helper function: checks if an atom is in any fully aromatic ring.
    def in_fully_aromatic_ring(atom):
        for ring in ringinfo.AtomRings():
            if atom.GetIdx() in ring:
                # If every atom in the ring is flagged aromatic, return True.
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    return True
        return False
    
    # Helper function: check if any immediate neighbor (except a given exclude) is aromatic.
    def neighbor_has_aromatic(atom, exclude_idx):
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == exclude_idx:
                continue
            if nbr.GetIsAromatic():
                return True
            if in_fully_aromatic_ring(nbr):
                return True
        return False

    # Examine each nitrile match.
    # Each match is a tuple: the first index is the nitrile carbon, the second is the nitrile nitrogen.
    for match in nitrile_matches:
        nitrile_carbon = mol.GetAtomWithIdx(match[0])
        nitrile_nitrogen = mol.GetAtomWithIdx(match[1])
        
        # Rule 1: The nitrile carbon must NOT be aromatic.
        if nitrile_carbon.GetIsAromatic():
            continue
        
        # Get the only substituent on the nitrile carbon (other than the nitrile nitrogen).
        # A nitrile carbon normally has degree 2, so we expect exactly one other neighbor.
        neighbors = [atom for atom in nitrile_carbon.GetNeighbors() if atom.GetIdx() != nitrile_nitrogen.GetIdx()]
        if not neighbors:
            continue  # not a typical nitrile
        
        substituent = neighbors[0]
        
        # Rule 2: The substituent should not be flagged as aromatic
        # either by its own aromatic flag or by being in a fully aromatic ring.
        if substituent.GetIsAromatic() or in_fully_aromatic_ring(substituent):
            continue
        
        # Rule 3: Check the immediate environment of the substituent.
        # Look at atoms directly attached to the substituent (except back to nitrile carbon);
        # if any such atom is aromatic we disqualify this nitrile group.
        if neighbor_has_aromatic(substituent, nitrile_carbon.GetIdx()):
            continue
        
        # Rule 4 (optional): Ensure that the bond between nitrile carbon and its substituent is a single bond.
        bond = mol.GetBondBetweenAtoms(nitrile_carbon.GetIdx(), substituent.GetIdx())
        if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
            continue

        # If we reach here for a nitrile match, we classify the nitrile group as aliphatic.
        return True, "Contains nitrile group with aliphatic environment"

    # If we have found nitrile groups but none satisfy our criteria, we return as false.
    return False, "Nitrile group(s) not found in a clean aliphatic environment"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        "CC(C)CC#N",                # isovaleronitrile, should pass.
        "CC#N",                     # acetonitrile, should pass.
        "ClC(Cl)C#N",               # Dichloroacetonitrile, should pass.
        "C\\C(=C\\CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)C#N",  # rhodiocyanoside A, borderline
        "O=C1CCC(C1)C#N",           # 3-oxocyclopentanecarbonitrile, should pass.
        "COC1=CC(=O)N(C)[C@@H](O)[C@@]1(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)C#N",  # acalyphin, should pass
        "c1ccccc1C#N",             # benzonitrile: aromatic nitrile, should fail.
    ]
    for smi in test_smiles:
        result, reason = is_aliphatic_nitrile(smi)
        print(f"SMILES: {smi}\n-> {result}, {reason}\n")