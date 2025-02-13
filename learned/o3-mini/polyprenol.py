"""
Classifies: CHEBI:26199 polyprenol
"""
#!/usr/bin/env python
"""
Classifies: polyprenol – any member of the class of prenols possessing the general formula
H-[CH2C(Me)=CHCH2]nOH in which n>=2 (i.e. more than one isoprene unit). 

This classifier requires that:
  • The molecule is acyclic (a linear isoprenoid chain).
  • The molecule is built exclusively from C, H and O.
  • It contains at least one hydroxyl (-OH) group, preferably at a terminal (primary) carbon.
  • The molecule contains a sufficient number of carbons (here, at least 10).
  • It displays at least two repeating isoprene-like units as detected by one or more SMARTS patterns.
  
We define three types of isoprene-like fragments:
    • Terminal (alpha) unit: [OH]CC([CH3])=C
    • Internal unit: [CH2]C([CH3])=C[CH2]
    • Terminal (omega) unit: C([CH3])=C[CH2]O
In our algorithm the sum of matches for these three patterns must be >=2.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    
    A polyprenol should:
       • Be acyclic (a linear chain, with no rings)
       • Be composed solely of C, H, and O (to avoid related phosphate, sulfur or acid compounds)
       • Have at least one hydroxyl group (-OH) that is located terminally (i.e. on a primary carbon)
       • Have a sufficiently long carbon chain (at least 10 C atoms)
       • Possess at least two isoprene–like repeating fragments.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a polyprenol, else False.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check that molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; expected a linear isoprenoid chain"
    
    # 2. Exclude molecules that contain atoms other than C, H, or O.
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains atom {atom.GetSymbol()} not in C, H, O"
    
    # 3. Check for at least one hydroxyl (-OH) group.
    hydroxyl_pat = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pat):
        return False, "No hydroxyl (-OH) group found"
    
    # 4. Look for a terminal hydroxyl group.
    terminal_oh_found = False
    for atom in mol.GetAtoms():
        # Check if atom is an -OH group (oxygen with one hydrogen)
        if atom.GetAtomicNum() == 8:
            # Count explicit hydrogens (or use GetTotalNumHs)
            if atom.GetTotalNumHs() >= 1:
                # Get attached heavy neighbors
                neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
                if len(neighbors) == 1:
                    neigh = neighbors[0]
                    # Ensure the neighbor is a carbon.
                    if neigh.GetAtomicNum() == 6:
                        # Count how many carbon neighbors that carbon has (excluding the O)
                        c_neighbors = [n for n in neigh.GetNeighbors() if n.GetAtomicNum() == 6]
                        # For a terminal hydroxyl the carbon should be primary (only one carbon neighbor).
                        if len(c_neighbors) == 1:
                            terminal_oh_found = True
                            break
    if not terminal_oh_found:
        return False, "No terminal hydroxyl (-OH) group on a primary carbon found"
    
    # 5. Check that there are a sufficient number of carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 10:
        return False, f"Too few carbons ({len(carbons)}) for a polyprenol structure"
    
    # (Optional) Check molecular weight; polyprenols usually have MW >150 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a polyprenol"
    
    # 6. Define SMARTS patterns for isoprene-like repeating fragments.
    # Terminal alpha unit: hydroxyl at beginning.
    pattern_alpha = Chem.MolFromSmarts("[OH]CC([CH3])=C")
    # Internal repeating unit.
    pattern_internal = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]")
    # Terminal omega unit: hydroxyl at end.
    pattern_omega = Chem.MolFromSmarts("C([CH3])=C[CH2]O")
    
    # Count matches for each pattern.
    matches_alpha = mol.GetSubstructMatches(pattern_alpha)
    matches_internal = mol.GetSubstructMatches(pattern_internal)
    matches_omega = mol.GetSubstructMatches(pattern_omega)
    total_units = len(matches_alpha) + len(matches_internal) + len(matches_omega)
    
    if total_units < 2:
        return False, f"Only {total_units} isoprene-like unit(s) detected (need at least 2)"
    
    return True, f"Polyprenol detected with {total_units} isoprene-like unit(s) and a terminal -OH group"

# Example usage:
if __name__ == "__main__":
    examples = {
        "(E,E,E)-geranylgeraniol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "Bionectin F": "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CO)/C)/C)C)C)C)C)C)C)C",
        "Dolichol-19": "OCC[C@H](CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC\\C=C(\\CC/C=C(/CCC=C(C)C)\\C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C",
        "Glisoprenin E": "OC(C(O)CC/C(=C/CO)/C)(CC/C=C(/CC/C=C(/CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC(O)C(O)(C)C)C)C)C)C)\\C)\\C)C",
        "geraniol": "CC(C)=CCC\\C(C)=C\\CO",
        "Dolichol-18": "OCC[C@H](CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC\\C=C(\\CC/C=C(/CCC=C(C)C)\\C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)",
        "Glisoprenin A": "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(=CCCC(=CCCC(=CCCC(=CCO)C)C)C)C)C)C)C)C",
        "(2-trans,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CO",
        "farnesol": "[H]C(CO)=C(C)CCC([H])=C(C)CCC=C(C)C",
        "SCH 66878": "OC(C(O)CC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC=C(C)C)C)C)C)C)C)C)C)C)C)C",
        "all-trans-octaprenol": "C(/C=C(/CC\\C=C(\\CC/C=C(\\C)/CCC=C(C)C)/C)\\C)C\\C(\\C)=C\\CC\\C(\\C)=C\\CC/C(/C)=C/CC\\C(=C\\CO)\\C",
        "ditrans,polycis-undecaprenol": "C(/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)/C)/C)/C)/C)/C)/C)/C)O",
        "Gymnoprenol A10": "OC(C(O)CO)(CC/C=C(/CC/C=C(/CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC=C(C)C)C)C)C)C)C)C)\\C)\\C)C",
        "(2-cis,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CO",
        "Glisoprenin D": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)/C)C)C)C)C",
        "all-trans-undecaprenol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "nerol": "C(=C\\CO)(\\CCC=C(C)C)/C",
        "(2-cis,6-cis)-farnesol": "CC(C)=CCC\\C(C)=C/CC\\C(C)=C/CO",
        "(2E,6E,10Z,14Z,18Z,22E)-3,7,11,15,19,23,27,31,35-Nonamethylhexatriaconta-2,6,10,14,18,22-hexaen-1-ol": "OC\\C=C(\\CC\\C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC\\C=C(\\CCCC(CCCC(CCCC(C)C)C)C)/C)/C)/C)/C)/C)/C",
        "geranylgeraniol": "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCO",
        "(2-trans,6-cis)-farnesol": "CC(C)=CCC\\C(C)=C/CC\\C(C)=C\\CO",
        "Glisoprenin F": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
        "all-trans-hexaprenol": "C(/C=C(/CCC=C(C)C)\\C)C\\C(\\C)=C\\CC\\C(\\C)=C\\CC/C(/C)=C/CC\\C(=C\\CO)\\C",
        "SCH 60057": "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
        "solanesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO"
    }
    
    for name, s in examples.items():
        result, explanation = is_polyprenol(s)
        print(f"{name}:\n  Result: {result}\n  Explanation: {explanation}\n")