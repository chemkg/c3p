"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid

A glycosphingolipid is defined as a carbohydrate-containing derivative of a sphingoid or ceramide.
The carbohydrate (typically in a pyranose sugar ring) is attached by a glycosidic linkage to O-1 of the sphingoid.
Heuristics used:
  - The molecule must be valid.
  - It must contain at least one pyranose sugar ring (a 6-membered ring with one oxygen).
  - It must contain an amide group (NC(=O)) typical of a ceramide.
  - At least one of the sugar ring carbons must have an exocyclic oxygen that connects further to a non-sugar atom,
    which we use as a crude proxy for a glycosidic linkage.
  - There should be enough non-ring aliphatic carbons to be consistent with a lipid tail.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule satisfies the criteria for a glycosphingolipid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Look for at least one pyranose ring.
    # The SMARTS pattern below represents a saturated six-membered ring (pyranose)
    # It requires one sp3 oxygen in the ring and five sp3 carbons.
    sugar_smarts = "[OX2r6]1[CX4r6][CX4r6][CX4r6][CX4r6][CX4r6]1"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No pyranose sugar ring found"
    
    # Heuristic 2: Look for an amide group typical of ceramides.
    # We need an "NC(=O)" fragment.
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group (typical of ceramide) found"
        
    # Heuristic 3: Look for a glycosidic linkage.
    # In glycosphingolipids the sugar is connected to the sphingoid through an oxygen that is bound
    # to a sugar-ring carbon (the anomeric carbon).  Many sugars have hydroxyl substituents;
    # here we look for a carbon atom in a detected sugar ring that has an oxygen neighbor which itself
    # is connected to an atom not in the sugar ring. This extra connectivity is used as a proxy for a glycosidic bond.
    glyco_linkage_found = False
    for match in sugar_matches:
        # The 'match' is a tuple of atom indices that make up the ring.
        sugar_indices = set(match)
        # Loop over carbon atoms in the ring (skip the ring oxygen)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            # Look at neighbors of this carbon; we want an oxygen that is exocyclic.
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == "O" and nbr.GetIdx() not in sugar_indices:
                    # To increase the chance that this oxygen is part of a glycosidic linkage,
                    # we check if this exocyclic oxygen is further connected to another atom that is also not in the sugar.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() != idx and nbr2.GetIdx() not in sugar_indices:
                            glyco_linkage_found = True
                            break
                    if glyco_linkage_found:
                        break
            if glyco_linkage_found:
                break
        if glyco_linkage_found:
            break
    if not glyco_linkage_found:
        return False, "No exocyclic oxygen (indicative of glycosidic linkage) detected from a sugar ring carbon"
    
    # Heuristic 4 (optional): Check for a sufficiently long aliphatic (lipid) chain.
    # Here we simply count non-ring carbon atoms.
    aliphatic_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if len(aliphatic_carbons) < 10:
        return False, "Not enough aliphatic carbons to indicate a lipid tail"
    
    return True, "Contains a pyranose sugar ring with an exocyclic oxygen (glycosidic linkage), an amide group, and a lipid tail"

# Example usage (can be removed or commented out when used as a module):
if __name__ == "__main__":
    # One example SMILES string: N-(2-hydroxyhexanoyl)-beta-D-galactosylsphinganine
    smiles_example = "CCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)C(O)CCCC"
    result, reason = is_glycosphingolipid(smiles_example)
    print("Result:", result)
    print("Reason:", reason)