"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid

A glycosphingolipid is a carbohydrate‐containing derivative of a sphingoid or ceramide.
The carbohydrate (typically in a pyranose sugar ring) is attached by a glycosidic linkage to O‑1 of the sphingoid.
Heuristics used in this improved version:
  - The molecule must be valid.
  - It must contain at least one pyranose sugar ring (six‐membered ring with one oxygen and five carbons).
  - It must have either a ceramide “amide” moiety (NC(=O)) or a sphingoid‐like portion (primary amine attached 
    to a stereogenic carbon bearing a –CH2OH group, detected via [NX3;H2][C@@H](CO) or [NX3;H2][C@H](CO)).
  - At least one sugar ring carbon must bear an exocyclic oxygen connected externally (a crude proxy for a glycosidic linkage).
  - There should be enough non‐ring aliphatic carbons (≥10) to indicate a lipid tail.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule meets our criteria for a glycosphingolipid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Look for at least one pyranose sugar ring.
    # SMARTS for a saturated six‐membered ring that contains one oxygen and five carbons.
    sugar_smarts = "[OX2r6]1[CX4r6][CX4r6][CX4r6][CX4r6][CX4r6]1"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No pyranose sugar ring found"
    
    # Heuristic 2: Look for either an amide group (as often seen in ceramides)
    # or a sphingoid base signature (with a primary amine attached to a stereogenic carbon with a CH2OH).
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    sphingoid_pattern1 = Chem.MolFromSmarts("[NX3;H2][C@@H](CO)")
    sphingoid_pattern2 = Chem.MolFromSmarts("[NX3;H2][C@H](CO)")
    if not (mol.HasSubstructMatch(amide_pattern) or 
            mol.HasSubstructMatch(sphingoid_pattern1) or 
            mol.HasSubstructMatch(sphingoid_pattern2)):
        return False, "Neither ceramide amide group nor sphingoid base signature (N[C](CO)) found"
    
    # Heuristic 3: Check for a glycosidic linkage.
    # We require that at least one sugar-ring carbon (excluding the ring oxygen)
    # has an exocyclic oxygen that connects to at least one other non-sugar atom.
    glyco_linkage_found = False
    for match in sugar_matches:
        # match is a tuple of atom indices that constitute the sugar ring.
        sugar_indices = set(match)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # Consider only carbons in the sugar (exclude the ring oxygen)
            if atom.GetSymbol() != "C":
                continue
            # For each neighbor of this carbon check if it is an oxygen not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == "O" and (nbr.GetIdx() not in sugar_indices):
                    # Check that this oxygen is connected further outside the sugar.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() != idx and (nbr2.GetIdx() not in sugar_indices):
                            glyco_linkage_found = True
                            break
                    if glyco_linkage_found:
                        break
            if glyco_linkage_found:
                break
        if glyco_linkage_found:
            break
    if not glyco_linkage_found:
        return False, "No exocyclic oxygen (indicative of a glycosidic linkage) detected from a sugar ring carbon"
    
    # Heuristic 4: Check for a sufficiently long aliphatic (lipid) chain.
    # Count all non-ring carbon atoms (atomic number 6) that are not in any ring.
    aliphatic_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if len(aliphatic_carbons) < 10:
        return False, "Not enough aliphatic carbons to indicate a lipid tail"
    
    return True, ("Contains a pyranose sugar ring with an exocyclic oxygen (glycosidic linkage), "
                  "and a sphingoid/ceramide moiety (amide or N[C](CO) pattern) along with a lipid tail")

# Example usage (can be removed when used as a module):
if __name__ == "__main__":
    # An example: N-(2-hydroxyhexanoyl)-beta-D-galactosylsphinganine
    smiles_example = "CCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)C(O)CCCC"
    result, reason = is_glycosphingolipid(smiles_example)
    print("Result:", result)
    print("Reason:", reason)