"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: Aliphatic Nitrile
Definition: Any nitrile derived from an aliphatic compound.
This program checks for the presence of a nitrile (C≡N) group and then
verifies that the nitrile group is attached to an aliphatic (carbon-only, sp³)
environment (i.e. not directly attached to an aromatic ring).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile must contain at least one nitrile group (a carbon triple-bonded to nitrogen)
    in which the nitrile carbon is not aromatic and is attached directly to an aliphatic carbon (sp³),
    not belonging to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for a nitrile group: an sp carbon triple bonded to a nitrogen.
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group found in the molecule"

    # Get ring information for later queries (to check aromatic rings)
    ringinfo = mol.GetRingInfo()

    # Iterate over all nitrile groups
    for match in nitrile_matches:
        # match is a tuple: first index is nitrile carbon, second is nitrile nitrogen
        nitrile_carbon = mol.GetAtomWithIdx(match[0])
        nitrile_nitrogen = mol.GetAtomWithIdx(match[1])

        # Rule 1: The nitrile carbon must not be aromatic.
        if nitrile_carbon.GetIsAromatic():
            # Skip this match if nitrile carbon is aromatic
            continue

        # Rule 2: Get the substituent attached to the nitrile carbon other than the nitrile N.
        neighbors = [atom for atom in nitrile_carbon.GetNeighbors() if atom.GetIdx() != nitrile_nitrogen.GetIdx()]
        if not neighbors:
            continue  # no substituent found; not a typical nitrile

        # For our purposes, we require that at least one substituent is a carbon with an aliphatic (sp³) environment.
        substituent = neighbors[0]
        # Check that the substituent is carbon
        if substituent.GetAtomicNum() != 6:
            continue
        # Check that the substituent is sp3-hybridized (a common trait of aliphatic chains)
        if substituent.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        # Check that the substituent is not part of any aromatic ring
        in_aromatic_ring = False
        for ring in ringinfo.AtomRings():
            if substituent.GetIdx() in ring:
                # If every atom in the ring is flagged as aromatic then consider it aromatic
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    in_aromatic_ring = True
                    break
        if in_aromatic_ring:
            continue

        # If one nitrile group passes all checks, we classify the molecule as an aliphatic nitrile.
        return True, "Contains nitrile group with aliphatic environment"

    # If no nitrile group passes the criteria then report failure.
    return False, "Nitrile group(s) not found in an aliphatic environment"

# Example usage (testing a few cases):
if __name__ == "__main__":
    test_smiles = [
        "CC#N",               # acetonitrile (aliphatic)
        "CC(C)CC#N",          # isovaleronitrile (aliphatic)
        "c1ccccc1C#N",        # benzonitrile (aromatic nitrile, should fail)
        "ClC(Cl)C#N",         # Dichloroacetonitrile (aliphatic)
        "COC1=CC(=O)N(C)[C@@H](O)[C@@]1(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)C#N"  # acalyphin
    ]
    for smi in test_smiles:
        result, reason = is_aliphatic_nitrile(smi)
        print(f"SMILES: {smi}\n-> {result}, {reason}\n")