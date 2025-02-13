"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: tannin
Definition: 'Any of a group of astringent polyphenolic vegetable principles or compounds, chiefly complex glucosides of catechol and pyrogallol.'
Heuristic: Require at least 2 aromatic rings and at least 2 occurrences of catechol- or pyrogallol-type substructures.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    A tannin is a polyphenolic compound with multiple aromatic rings decorated by hydroxyl groups,
    such as catechol (dihydroxybenzene) or pyrogallol (trihydroxybenzene) moieties. Additionally, tannins 
    often have glycosidic components.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a tannin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count aromatic rings: use ring information and check if all atoms in a ring are aromatic.
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    if aromatic_ring_count < 2:
        return False, f"Only {aromatic_ring_count} aromatic ring(s) found; tannins typically have several polyphenolic rings"
    
    # Define SMARTS patterns for catechol and pyrogallol motifs.
    # Note: These patterns are simplified representations.
    catechol_smarts = "c1ccc(O)c(O)c1"      # benzene with two hydroxyl groups in ortho arrangement
    pyrogallol_smarts = "c1cc(O)c(O)c(O)c1"  # benzene with three hydroxyl groups
    catechol_pattern = Chem.MolFromSmarts(catechol_smarts)
    pyrogallol_pattern = Chem.MolFromSmarts(pyrogallol_smarts)
    
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pattern)
    total_polyphenol_units = len(catechol_matches) + len(pyrogallol_matches)
    
    if total_polyphenol_units < 2:
        return False, f"Found only {total_polyphenol_units} catechol/pyrogallol unit(s); tannins typically contain multiple such polyphenolic motifs"
    
    # Optional criterion: check molecular weight. Tannins are not very small, often >300 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is lower than expected for a tannin"
    
    # If all heuristic criteria are met, we classify it as a tannin.
    return True, f"Contains {aromatic_ring_count} aromatic ring(s) and {total_polyphenol_units} polyphenolic unit(s) with a molecular weight of {mol_wt:.1f} Da, consistent with tannins"

# Example usage:
if __name__ == "__main__":
    test_smiles = "O=C(O)C1=CC(OC)=C(O[C@@H]2O[C@@H]([C@@H](OC)[C@@H]([C@H]2O)O)CO)C(=C1)CC=C(C)C"  # Conoideoglucoside C
    result, reason = is_tannin(test_smiles)
    print(result, reason)