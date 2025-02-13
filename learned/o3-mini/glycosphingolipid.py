"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid
A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide.
The carbohydrate residue is attached by a glycosidic linkage to O-1 of the sphingoid.
This routine uses heuristic SMARTS patterns to look for a sugar ring, a ceramide-like (amide) motif,
and evidence for a glycosidic bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    Heuristics used:
      - The molecule must be a valid structure.
      - It must contain at least one pyranose sugar ring.
      - It must contain an amide group (NC(=O)) typical of ceramides.
      - At least one oxygen within the sugar ring must be linked to an atom outside the ring,
        as expected for a glycosidic linkage from the sugar to the sphingoid.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if it fulfills our criteria for a glycosphingolipid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Heuristic 1: Look for a pyranose ring.
    # This SMARTS roughly matches a pyranose ring with defined hydroxyl groups.
    # (Stereochemistry is specified in many glycosidic examples, but we ignore some variation.)
    sugar_smarts = "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No pyranose sugar ring found"
    
    # Heuristic 2: Look for a ceramide-like (amide) group.
    # Ceramides typically contain an amide bond (NC(=O)).
    amide_pattern = Chem.MolFromSmarts("N[C;!$(C(=O))]C(=O)")  # a simple amide search
    # In practice, many ceramides will show an "N-C(=O)" moiety; we relax the pattern slightly.
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("NC(=O)")):
        return False, "No amide group (typical of ceramides) found"
    
    # Optionally, check for a long carbon chain as a crude test for a lipid tail.
    # We count the number of aliphatic carbons.
    aliphatic_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.IsInRing() == False]
    if len(aliphatic_carbons) < 10:
        return False, "Not enough aliphatic carbons to indicate a lipid tail"
    
    # Heuristic 3: Check that at least one oxygen in the sugar ring is linked (via a single bond)
    # to an atom that is not part of the sugar. This is a crude check for a glycosidic linkage.
    linkage_found = False
    for match in sugar_matches:
        # get set of atom indices in the sugar ring match
        sugar_indices = set(match)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "O":
                # Check neighbors of this oxygen; if any neighbor is outside the sugar ring, it could be the glycosidic linkage.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in sugar_indices:
                        linkage_found = True
                        break
            if linkage_found:
                break
        if linkage_found:
            break
    if not linkage_found:
        return False, "No oxygen outside the sugar ring detected as a glycosidic linkage"

    # Passed all checks: classify as glycosphingolipid.
    return True, "Contains a pyranose sugar ring attached via glycosidic linkage to a ceramide-like (sphingoid) portion"

# Example usage (you can remove or comment out these lines if using as a library):
if __name__ == "__main__":
    # Test with one of the provided examples:
    smiles_example = "CCCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)C(O)CCCC"
    result, reason = is_glycosphingolipid(smiles_example)
    print("Result:", result)
    print("Reason:", reason)