"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: Unsaturated Fatty Acid
Definition: An unsaturated fatty acid is defined here as an acyclic molecule that:
  1. Contains exactly one terminal free carboxylic acid group (i.e. –C(=O)[OH]). The carboxyl carbon
     must be acyclic and only bonded to exactly one carbon (which is non‐aromatic).
  2. Contains only C, H, and O atoms.
  3. Contains at least one C=C or C#C bond (unsaturation) outside of the carboxyl group.
  4. Has a minimum chain size (at least 5 carbon atoms).
By enforcing that the acid group is a free –COOH (and not part of an ester or other motifs),
and by limiting the allowed elements, many false positives (e.g. phospholipids or acetates) are avoided.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    Criteria:
      1. The molecule must be acyclic.
      2. The molecule must contain only C, H, and O atoms.
      3. The molecule must have exactly one terminal free carboxylic acid group.
         (This is defined by a [CX3](=O)[OH] pattern on a carbon that is not in a ring and that
         has exactly one carbon neighbor (which must be non-aromatic)).
      4. The molecule must contain at least one carbon–carbon double or triple bond (unsaturation)
         outside of the carboxyl group.
      5. As a simple size check, the overall number of carbons should be at least 5.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as an unsaturated fatty acid; False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that contain rings – fatty acids are defined here as acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s); fatty acids are typically acyclic"
    
    # Reject molecules containing atoms other than C, H, and O.
    allowed = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed:
            return False, "Molecule contains atoms other than C, H, and O; not a typical fatty acid"
    
    # Look for terminal free carboxylic acid groups.
    # Use SMARTS that strictly requires a free carboxyl: –C(=O)[OH]
    acid_smarts = "[CX3](=O)[OH]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    terminal_acid_matches = []
    acid_bond_idxs = set()  # store indices of bonds that are part of the acid group
    for match in acid_matches:
        # match[0] is the carboxyl carbon for the acid group
        acid_c = mol.GetAtomWithIdx(match[0])
        if acid_c.IsInRing():
            continue
        # Count carbon neighbors of the acid carbon (ignoring oxygens)
        carbon_neighbors = [nb for nb in acid_c.GetNeighbors() if nb.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1 and not carbon_neighbors[0].GetIsAromatic():
            terminal_acid_matches.append(match)
            # Record the bonds that form the carboxyl group so we can exclude them
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                for bond in atom.GetBonds():
                    # If both bond endpoints are in the current acid match, record it.
                    if bond.GetBeginAtomIdx() in match and bond.GetEndAtomIdx() in match:
                        acid_bond_idxs.add(bond.GetIdx())
    
    if len(terminal_acid_matches) != 1:
        return False, ("Molecule does not have exactly one terminal free carboxylic acid group "
                       f"(found {len(terminal_acid_matches)}); not classified as a fatty acid")
    
    # Check overall carbon count.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbons to be a fatty acid"
    
    # Look for unsaturation (C=C or C#C bonds) outside of the carboxyl group.
    has_unsaturation = False
    unsat_types = set()
    for bond in mol.GetBonds():
        # Skip bonds that are part of the terminal acid group
        if bond.GetIdx() in acid_bond_idxs:
            continue
        if bond.GetBondType() == Chem.BondType.DOUBLE or bond.GetBondType() == Chem.BondType.TRIPLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # We only count if both atoms are carbon.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                has_unsaturation = True
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    unsat_types.add("double")
                else:
                    unsat_types.add("triple")
    
    if not has_unsaturation:
        return False, ("Molecule appears to have a free carboxylic acid group but does not contain "
                       "any carbon–carbon double or triple bonds (unsaturation)")
    
    # Prepare a description of the unsaturation found.
    if "double" in unsat_types and "triple" in unsat_types:
        unsat_detail = "carbon–carbon double and triple bond(s)"
    elif "double" in unsat_types:
        unsat_detail = "carbon–carbon double bond(s)"
    else:
        unsat_detail = "carbon–carbon triple bond(s)"
        
    return True, f"Molecule is a fatty acid and contains unsaturation: {unsat_detail}"

# Example usage:
if __name__ == "__main__":
    # Test with trans-vaccenic acid, one of the provided examples.
    test_smiles = "CCCCCC\\C=C\\CCCCCCCCCC(O)=O"
    result, reason = is_unsaturated_fatty_acid(test_smiles)
    print(result, reason)