"""
Classifies: CHEBI:39437 tocol
"""
#!/usr/bin/env python3
"""
Classifies: tocol
A tocol is defined as "A chromanol with a chroman-6-ol skeleton that is substituted at position 2 by a 
saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units".

We use a heuristic approach:
  (1) Search for a chromanol core via a SMARTS pattern that is flexible enough to match variants.
  (2) Among the atoms in the core, look for one (typically sp3 carbon) with a substituent (side chain) 
      that is not part of the core.
  (3) Extract the side chain (via a DFS where we do not traverse core atoms), and check that it is:
       - Acyclic (no rings).
       - Composed entirely of carbon atoms.
       - Contains between ~12 and 17 carbon atoms (consistent with three isoprene units).
       - Has either 0 (saturated) or 3 (triply unsaturated) double bonds.
If these checks pass, the molecule is classified as a tocol.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines whether the molecule specified by the SMILES string is a tocol.
    
    A tocol is defined as a chromanol (chroman-6-ol) skeleton substituted at position 2 by an
    appropriate hydrocarbon (isoprenoid) chain.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): A tuple: (True, reason) if tocol; otherwise (False, reason)
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more flexible SMARTS for the chromanol (chroman-6-ol) core.
    # The pattern here finds a fused bicyclic system: a saturated 6-membered ring with an oxygen
    # (the dihydropyran) fused to an aromatic ring which must have at least one hydroxyl (OH).
    # The pattern “O1CCc2cc(O)ccc2C1” should match many vitamin E core variants.
    chromanol_smarts = "O1CCc2cc(O)ccc2C1"
    chromanol_query = Chem.MolFromSmarts(chromanol_smarts)
    if chromanol_query is None:
        return False, "Error in chromanol SMARTS pattern"
    
    # Search for the chromanol substructure in the molecule
    core_match = mol.GetSubstructMatch(chromanol_query)
    if not core_match:
        return False, "Chromanol (chroman-6-ol) core not found"
    
    # Store the indices of the core atoms from the match.
    core_atoms = set(core_match)
    
    # Identify candidate attachment points within the core.
    # We expect the side chain to be attached to one of the carbon atoms of the saturated ring.
    candidate_indices = []
    for idx in core_atoms:
        atom = mol.GetAtomWithIdx(idx)
        # restrict to carbon atoms that are not aromatic:
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            # Look for a neighbor that is not in the core and is a carbon (the start of a side chain)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() not in core_atoms and nb.GetAtomicNum() == 6:
                    candidate_indices.append(idx)
                    break

    if not candidate_indices:
        return False, "Side chain attachment point not found on chromanol core"
    
    # For each candidate attachment point, try to extract a side chain subgraph.
    # We use DFS starting from the neighbor (outside of core) attached to the candidate.
    def extract_side_chain(start_idx):
        # DFS: do not cross back into core and only allow carbon atoms.
        visited = set()
        stack = [start_idx]
        side_chain = set()
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            side_chain.add(curr)
            curr_atom = mol.GetAtomWithIdx(curr)
            # If any side-chain atom is part of a ring, flag an error (acyclic expected)
            if curr_atom.IsInRing():
                return None, "Side chain is cyclic, expected a linear hydrocarbon chain"
            for nb in curr_atom.GetNeighbors():
                if nb.GetIdx() in core_atoms or nb.GetIdx() in visited:
                    continue
                # Ensure the neighbor is a carbon; if not, error.
                if nb.GetAtomicNum() != 6:
                    return None, "Side chain contains non-carbon atoms"
                stack.append(nb.GetIdx())
        return side_chain, None
    
    valid_side_chain_found = False
    side_reason = ""
    for candidate in candidate_indices:
        cand_atom = mol.GetAtomWithIdx(candidate)
        # Find neighbor(s) of the candidate outside the core to start the side chain.
        side_chain_found = False
        for nb in cand_atom.GetNeighbors():
            if nb.GetIdx() in core_atoms or nb.GetAtomicNum() != 6:
                continue
            # Attempt to extract the side chain starting from this neighbor.
            side_chain, err = extract_side_chain(nb.GetIdx())
            if err:
                side_chain_found = False
                side_reason = err
                continue
            # Count the number of carbon atoms in the chain.
            chain_carbons = sum(1 for idx in side_chain if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # Check the expected carbon count for an isoprenoid chain of three units (approximately 12-17 carbons).
            if chain_carbons < 12 or chain_carbons > 17:
                side_reason = f"Side chain has {chain_carbons} carbons; expected between 12 and 17 for isoprenoid units"
                continue
            # Count the number of double bonds in the side chain.
            chain_double_bonds = 0
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if a1 in side_chain and a2 in side_chain:
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        chain_double_bonds += 1
            if chain_double_bonds not in (0, 3):
                side_reason = f"Side chain has {chain_double_bonds} double bonds; expected 0 (saturated) or 3 (triply unsaturated)"
                continue
            # If we pass all tests, we have found a valid side chain.
            valid_side_chain_found = True
            side_chain_found = True
            break  # no need to check further neighbors for this candidate
        if valid_side_chain_found:
            break

    if not valid_side_chain_found:
        # Prioritize a meaningful message if available; else provide a generic failure message.
        if side_reason:
            return False, side_reason
        return False, "No appropriate side chain found attached to the chromanol core"
    
    return True, "Molecule contains a chromanol core with an appropriate isoprenoid side chain (tocol)"
    
# For testing when run as a script:
if __name__ == "__main__":
    # Example: beta-tocopherol
    test_smiles = "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@]1(C)CCc2c(C)c(O)cc(C)c2O1"
    result, reason = is_tocol(test_smiles)
    print(f"Test SMILES: {test_smiles}")
    print(f"Is tocol? {result}. Reason: {reason}")