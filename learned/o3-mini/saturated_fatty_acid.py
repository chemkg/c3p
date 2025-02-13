"""
Classifies: CHEBI:26607 saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Saturated fatty acid

A saturated fatty acid is defined as a fatty acid containing exactly one carboxylic acid group,
an open-chain aliphatic moiety with no carbon–carbon multiple bonds (the only allowed C=O is in the acid group)
and only C, H, and O atoms. We further disallow additional carbonyl (keto) groups in a long-chain (>=7 C)
and also reject cases where an -OH group is present directly on the alpha carbon.
Known to produce adverse biological effects when ingested to excess.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    The molecule must:
      - Contain only C, H, and O atoms.
      - Have exactly one carboxylic acid group (defined as [CX3](=O)[OX2H]).
      - Contain no rings.
      - Have no carbon–carbon double or triple bonds.
      - In molecules with seven or more carbons, not display additional carbonyl (C=O) groups 
        beyond that in the carboxylic acid.
      - Not have an -OH substituent on the alpha carbon (the carbon directly bonded to the acid carbon).
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a saturated fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check allowed atoms: only C (6), H (1), O (8)
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Contains atoms other than C, H, and O"
    
    # Check that the molecule contains no rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structures"
    
    # Define the carboxylic acid SMARTS pattern and ensure exactly one match.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    if len(acid_matches) != 1:
        return False, "More than one carboxylic acid group found"
    
    # Save the indices of the acid group atoms.
    # The SMARTS is defined so that: acid_match[0][0] = acid carbon, acid_match[0][1] = carbonyl oxygen, acid_match[0][2] = hydroxyl oxygen.
    acid_match = acid_matches[0]
    acid_carbon_idx = acid_match[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    acid_carbony_ox_idx = acid_match[1]  # allowed carbonyl oxygen for the acid group

    # Check for carbon-carbon multiple bonds (double or triple) anywhere in the molecule.
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                return False, "Contains carbon-to-carbon multiple bonds"
    
    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for additional carbonyl groups (C=O) not part of the acid group.
    # We allow extra C=O only in molecules with fewer than 7 carbons.
    if total_carbons >= 7:
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                # look for a bond between a carbon and an oxygen
                ats = {a1.GetAtomicNum(), a2.GetAtomicNum()}
                if 6 in ats and 8 in ats:
                    # Check if this bond is the allowed acid carbonyl.
                    # Allowed if one end is the acid carbon and the other is the designated carbonyl oxygen.
                    idxs = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
                    if acid_carbon_idx in idxs and acid_carbony_ox_idx in idxs:
                        continue  # allowed acid carbonyl
                    else:
                        return False, "Contains additional carbonyl group(s) beyond the carboxylic acid"
    
    # Check that the acid group is terminal. Reject if the atom alpha to the acid carbon has an -OH substituent.
    # Identify the neighbor(s) of the acid carbon excluding the acid carbonyl oxygen.
    for neighbor in acid_carbon.GetNeighbors():
        if neighbor.GetIdx() == acid_carbony_ox_idx:
            continue  # skip the acid carbonyl oxygen
        # Check if this neighbor (the alpha carbon) bears a hydroxyl substituent.
        for subnbr in neighbor.GetNeighbors():
            # Looking for an oxygen attached via a single bond and with at least one hydrogen.
            if subnbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), subnbr.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if the O has any attached hydrogens.
                    if subnbr.GetTotalNumHs() > 0:
                        return False, "Contains an alpha-hydroxyl group adjacent to the carboxylic acid"
    
    return True, ("Molecule has a single carboxylic acid group, only C/H/O atoms, no rings, "
                  "no carbon–carbon unsaturation, and no disallowed additional carbonyl or alpha‐OH groups")

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "OCCCCCCCCCCCCCCCCCCCC(O)=O",  # 20-hydroxyicosanoic acid; should be True
        "CC[C@@H](O)C(O)=O",           # (R)-2-hydroxybutyric acid; should be False (alpha-OH)
        "CC(=O)CCC(O)=O",              # 4-oxopentanoic acid; small keto acid, allowed
        "O=C(CCCCCCCCC(O)=O)CCCC",      # 10-keto myristic acid; should be False (extra carbonyl in long chain)
    ]
    for sm in test_smiles:
        result, reason = is_saturated_fatty_acid(sm)
        print(f"SMILES: {sm} -> {result}, Reason: {reason}")