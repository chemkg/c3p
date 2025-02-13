"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: Aliphatic Nitrile
Definition: Any nitrile derived from an aliphatic compound.

This improved program checks for the presence of a nitrile (C≡N) group and then verifies that:
 1. The nitrile carbon is non‐aromatic.
 2. Its only substituent (the atom attached to the nitrile carbon besides nitrogen) is non‐aromatic
    and is sp³-hybridized.
 3. None of the immediate neighbors of that substituent (except for the nitrile carbon) display aromaticity.
These additional checks help to prevent false positives in cases where the nitrile group is
attached to a carbon that is part of unsaturated or aromatic motifs.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    We require that at least one nitrile group (C≡N) is present where the nitrile carbon
    is non‐aromatic and is bonded to a substituent that is aliphatic (non‐aromatic and sp3-hybridized)
    and its immediate neighbors (except the bond to the nitrile carbon) are also free from aromaticity.

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
    
    # SMARTS pattern for a nitrile group: a carbon triple bonded to a nitrogen
    nitrile_pattern = Chem.MolFromSmarts("[C]#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group found in the molecule"
    
    ringinfo = mol.GetRingInfo()  # to check ring aromaticity

    # Helper function: checks if an atom is in any fully aromatic ring.
    def in_fully_aromatic_ring(atom):
        for ring in ringinfo.AtomRings():
            if atom.GetIdx() in ring:
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    return True
        return False
    
    # Helper function: check if any immediate neighbor (besides a given exclude) is aromatic
    def neighbor_has_aromatic(atom, exclude_idx):
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == exclude_idx:
                continue
            if nbr.GetIsAromatic():
                return True
            if in_fully_aromatic_ring(nbr):
                return True
        return False

    # Iterate over each nitrile match.
    # Each match is a tuple: (nitrile carbon index, nitrile nitrogen index).
    for match in nitrile_matches:
        nitrile_carbon = mol.GetAtomWithIdx(match[0])
        nitrile_nitrogen = mol.GetAtomWithIdx(match[1])
        
        # Rule 1: The nitrile carbon must not be aromatic.
        if nitrile_carbon.GetIsAromatic():
            continue
        
        # Get the only substituent on the nitrile carbon (other than the nitrile nitrogen).
        neighbors = [atom for atom in nitrile_carbon.GetNeighbors() if atom.GetIdx() != nitrile_nitrogen.GetIdx()]
        if not neighbors:
            continue  # atypical nitrile
        
        substituent = neighbors[0]
        
        # Rule 2: The substituent must not be flagged as aromatic,
        # and it should not be in a fully aromatic ring.
        if substituent.GetIsAromatic() or in_fully_aromatic_ring(substituent):
            continue
        
        # Rule 2b: The substituent should be sp3 hybridized (i.e. aliphatic)
        if substituent.GetHybridization() != rdchem.HybridizationType.SP3:
            continue

        # Rule 3: Check the immediate neighbourhood of the substituent (except back to nitrile carbon).
        if neighbor_has_aromatic(substituent, nitrile_carbon.GetIdx()):
            continue
        
        # Rule 4: Ensure that the bond between the nitrile carbon and its substituent is a single bond.
        bond = mol.GetBondBetweenAtoms(nitrile_carbon.GetIdx(), substituent.GetIdx())
        if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
            continue

        return True, "Contains nitrile group with aliphatic environment"

    # If no nitrile group satisfies our stricter criteria then we return False.
    return False, "Nitrile group(s) not found in a clean aliphatic environment"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        "CC(C)CC#N",                # isovaleronitrile, should pass.
        "CC#N",                     # acetonitrile, should pass.
        "ClC(Cl)C#N",               # Dichloroacetonitrile, should pass.
        "C\\C(=C\\CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)C#N",  # rhodiocyanoside A, borderline
        "O=C1CCC(C1)C#N",           # 3-oxocyclopentanecarbonitrile, should pass.
        "CCCC(N)C#N",              # 2-aminopentanenitrile, should pass.
        "COC1=CC(=O)N(C)[C@@H](O)[C@@]1(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)C#N",  # acalyphin, should pass
        "c1ccccc1C#N",             # benzonitrile (aromatic nitrile), should fail.
        # False positives (should fail classification)
        "COCC(=O)N1[C@@H]([C@@H]([C@H]1C#N)C2=CC=C(C=C2)C#CCC3CCCC3)CO",
        "CCOC(=O)C1=CC(=NC2=C1C(=NN2CCC#N)C)C3=C(C=C(C=C3)OC)F",
    ]
    for smi in test_smiles:
        result, reason = is_aliphatic_nitrile(smi)
        print(f"SMILES: {smi}\n-> {result}, {reason}\n")