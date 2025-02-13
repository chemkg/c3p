"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: Unsaturated Fatty Acid
Definition: A fatty acid is defined here as an acyclic molecule with exactly one terminal (acyclic, non-aromatic) carboxyl group,
and an aliphatic carbon chain that contains at least one carbon–carbon double or triple bond (apart from the carbonyl).
The extra criteria (acyclic structure, minimum number of carbons, and non-aromatic attachment of the acid group)
help to avoid false positives.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    Criteria:
      1. The molecule must be acyclic.
      2. The molecule must have exactly one terminal (acyclic, non-aromatic) carboxylic acid group.
         (This is defined by looking for a C(=O)O or C(=O)[O-] group where the acid C is bonded to exactly one carbon neighbor that is not aromatic.)
      3. The molecule must contain at least one carbon–carbon double or triple bond outside of the carboxyl group.
      4. As a simple size check, the overall number of carbons should be at least 5.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as an unsaturated fatty acid, else False.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Fatty acids are aliphatic; reject molecules containing any rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s); fatty acids are typically acyclic"
    
    # Look for carboxyl groups (both neutral and deprotonated forms)
    acid_smarts = "[$(C(=O)O),$(C(=O)[O-])]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Filter for terminal acid groups:
    # The acid carbon (first atom in the match) must be acyclic and have exactly one carbon neighbor that is non-aromatic.
    terminal_acid_matches = []
    for match in acid_matches:
        acid_c = mol.GetAtomWithIdx(match[0])
        if acid_c.IsInRing():
            continue
        # Consider only carbon neighbors (ignore oxygens)
        carbon_neighbors = [nb for nb in acid_c.GetNeighbors() if nb.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1 and not carbon_neighbors[0].GetIsAromatic():
            terminal_acid_matches.append(match)
    
    if len(terminal_acid_matches) != 1:
        return False, ("Molecule does not have exactly one terminal carboxylic acid group "
                       f"(found {len(terminal_acid_matches)}), so not a fatty acid")
    
    # Check that the molecule has at least a minimum number of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbons to be a fatty acid"
    
    # To ensure that we are not counting the carbonyl in the acid group as unsaturation,
    # exclude bonds that are part of the terminal carboxyl group.
    acid_bond_idxs = set()
    for match in terminal_acid_matches:
        acid_idx = match[0]
        for bond in mol.GetAtomWithIdx(acid_idx).GetBonds():
            acid_bond_idxs.add(bond.GetIdx())
    
    # Now check for unsaturation: at least one C=C or C#C bond outside of those in the carboxyl group.
    has_unsaturation = False
    unsat_bonds = []  # collect unsaturation bonds for reporting
    for bond in mol.GetBonds():
        # Skip bonds that are part of the terminal acid group
        if bond.GetIdx() in acid_bond_idxs:
            continue
        # Check if the bond is a double or triple bond and that both atoms are carbon.
        if (bond.GetBondType() == Chem.BondType.DOUBLE or bond.GetBondType() == Chem.BondType.TRIPLE):
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                has_unsaturation = True
                unsat_bonds.append(bond)
    
    if not has_unsaturation:
        return False, ("Molecule appears to be a fatty acid (has a terminal carboxyl group) but does not contain "
                       "any C=C or C#C bonds (unsaturation)")
    
    # Prepare a message about the type of unsaturation found.
    if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in unsat_bonds) and any(bond.GetBondType() == Chem.BondType.TRIPLE for bond in unsat_bonds):
        unsat_detail = "carbon-carbon double and triple bond(s)"
    elif any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in unsat_bonds):
        unsat_detail = "carbon-carbon double bond(s)"
    else:
        unsat_detail = "carbon-carbon triple bond(s)"
    
    return True, f"Molecule is a fatty acid and contains unsaturation: {unsat_detail}"

# Example usage:
if __name__ == "__main__":
    # Test with trans-vaccenic acid, one of the provided examples.
    test_smiles = "CCCCCC\\C=C\\CCCCCCCCCC(O)=O"
    result, reason = is_unsaturated_fatty_acid(test_smiles)
    print(result, reason)