"""
Classifies: CHEBI:39437 tocol
"""
#!/usr/bin/env python3
"""
Classifies: tocol
A tocol is defined as "A chromanol with a chroman-6-ol skeleton that is substituted at position 2
by a saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units".
This program uses a heuristic substructure search:
  (1) It looks for a chromanol (chroman-6-ol) core.
  (2) It finds the side chain attached at the position 2 carbon.
  (3) It verifies that the side chain (ignoring the core atoms) is composed entirely of carbons,
      is acyclic, has a carbon count in an acceptable range (~12 to 17), and that it is either fully saturated
      (0 double bonds) or triply unsaturated (3 double bonds).
If these criteria are met, we classify the molecule as a tocol.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines whether the molecule specified by the SMILES string is a tocol.
    
    A tocol is defined as a chromanol (chroman-6-ol) skeleton substituted at position 2 by a 
    saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): A tuple where the first element is True if the molecule is classified
                     as a tocol, and False otherwise. The second element provides a reason.
    """
    # Parse SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the chromanol core.
    # The pattern is defined as a benzopyran ring (a six-membered aromatic ring fused to a saturated ring containing an oxygen)
    # bearing an OH group on the aromatic ring.
    # We label the first saturated carbon (immediately attached to the ring oxygen) as [C:2] because in vitamin E
    # derivatives that is the carbon (position 2) that carries the hydrocarbon side chain.
    #
    # Note: The pattern is kept somewhat simplified and may not cover every valid stereochemical or methyl-substituted variant.
    chromanol_smarts = "O1[C:2]CCc2ccc(O)cc2C1"
    chromanol_query = Chem.MolFromSmarts(chromanol_smarts)
    if chromanol_query is None:
        return False, "Error in chromanol SMARTS pattern"

    # Look for the chromanol substructure in the molecule.
    core_match = mol.GetSubstructMatch(chromanol_query)
    if not core_match:
        return False, "Chromanol (chroman-6-ol) core not found"

    # Our SMARTS pattern places the [C:2] as the second atom in the query.
    # In the match tuple, index 1 corresponds to that atom.
    pos2_idx = core_match[1]
    pos2_atom = mol.GetAtomWithIdx(pos2_idx)
    
    # Identify the side chain: the substituent(s) on the pos2_atom that are not part of the core match.
    core_atoms = set(core_match)
    side_chain_start = None
    for nb in pos2_atom.GetNeighbors():
        if nb.GetIdx() not in core_atoms and nb.GetAtomicNum() == 6:  # ensure it is a carbon
            side_chain_start = nb
            break
    if side_chain_start is None:
        return False, "Side chain not found at position 2"
    
    # Now we extract the subgraph (side chain) attached at the pos2 atom.
    # We do a simple depth-first search starting from the side_chain_start atom, excluding atoms that are known to be in the core.
    visited = set()
    to_visit = [side_chain_start.GetIdx()]
    side_chain_atoms = set()
    while to_visit:
        curr_idx = to_visit.pop()
        if curr_idx in visited:
            continue
        visited.add(curr_idx)
        side_chain_atoms.add(curr_idx)
        curr_atom = mol.GetAtomWithIdx(curr_idx)
        # We only continue traversal if the atom is not part of any ring (to ensure it is a simple chain)
        # Because by definition the side chain is a hydrocarbon chain.
        if curr_atom.IsInRing():
            # If a ring is encountered in the side chain, then the substituent is not a simple hydrocarbon chain.
            return False, "Side chain is cyclic, expected a linear hydrocarbon chain"
        # Traverse neighbors not in the core.
        for nb in curr_atom.GetNeighbors():
            if nb.GetIdx() not in core_atoms and nb.GetIdx() not in visited:
                # We only follow carbon atoms in the chain.
                if nb.GetAtomicNum() == 6:
                    to_visit.append(nb.GetIdx())
                else:
                    # If any heteroatom is found in the chain, it does not match the expected hydrocarbon chain.
                    return False, "Side chain contains non-carbon atoms"
    
    # Count the number of carbons in the side chain.
    chain_carbons = sum(1 for idx in side_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # For an isoprenoid chain of three isoprene units, one would expect roughly 15 carbons.
    # Allow a bit of margin (e.g., between 12 and 17 carbons).
    if chain_carbons < 12 or chain_carbons > 17:
        return False, f"Side chain has {chain_carbons} carbons; expected between 12 and 17 for three isoprenoid units"
    
    # Count the number of double bonds in the side chain.
    chain_double_bonds = 0
    # We iterate over all bonds in the molecule; if both atoms belong to the side chain subgraph then count.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in side_chain_atoms and a2 in side_chain_atoms:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                chain_double_bonds += 1
    
    # The definition permits either a saturated chain (0 double bonds) or a triply unsaturated chain (3 double bonds).
    if chain_double_bonds not in (0, 3):
        return False, f"Side chain has {chain_double_bonds} double bonds; expected 0 (saturated) or 3 (triply-unsaturated)"
    
    return True, "Molecule contains chromanol core with an appropriate isoprenoid side chain (tocol)"
    
# For testing purposes (if run as a script), one could include:
if __name__ == "__main__":
    # Example: beta-tocopherol
    test_smiles = "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@]1(C)CCc2c(C)c(O)cc(C)c2O1"
    result, reason = is_tocol(test_smiles)
    print(f"Test SMILES: {test_smiles}")
    print(f"Is tocol? {result}. Reason: {reason}")