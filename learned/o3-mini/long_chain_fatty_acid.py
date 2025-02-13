"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: Long-chain fatty acid (chain length ranging from C13 to C22)

A long-chain fatty acid is defined as a fatty acid with a chain length ranging 
from C13 to C22. In our implementation we require that the molecule contains a terminal 
carboxyl group, which we detect by finding a carbon atom (acid carbon) that:
  - Is bonded to exactly two oxygens: one via a double bond (carbonyl) and one via a single bond (hydroxyl or anionic oxygen)
  - Is attached to exactly one other carbon (making it “terminal”)
Then we “walk” the connected carbon chain (ignoring non-carbon atoms) while enforcing that:
  - The chain is acyclic (none of the chain carbons is in a ring)
  - The chain is unbranched (each chain carbon along the path extends in only one direction)
If the total number of carbons in this chain (including the carboxyl carbon) is between 13 and 22, 
the molecule is classified as a long-chain fatty acid.
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a long-chain fatty acid based on its SMILES string.

    A long-chain fatty acid is defined as having a terminal carboxyl group (C(=O)O or C(=O)[O-])
    whose acid carbon is attached to exactly one other carbon, and the longest linear carbon chain 
    (starting at that neighbor and including the acid carbon) has a length between 13 and 22. Additionally,
    the carbon chain must be acyclic and unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule qualifies as a long-chain fatty acid, False otherwise.
        str: Reason detailing the classification outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, iterate over atoms to find a candidate acid carbon
    # Our criteria for an acid (carboxyl) carbon are:
    #  1. It is a carbon (atomic num 6)
    #  2. It is bonded to exactly two oxygens: one by a double bond and one by a single bond.
    #  3. Aside from these oxygens, it is bonded to exactly one other carbon.
    candidate_found = False
    candidate_acid_atom = None
    candidate_chain_neighbor = None

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Find oxygen neighbors (ignoring any other atoms)
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # record the bond type as well
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                oxy_neighbors.append((nbr, bond.GetBondType()))
        if len(oxy_neighbors) != 2:
            continue
        
        # Check that one oxygen is double bonded and one is single bonded.
        dbl_bond = False
        sgl_bond = False
        for nbr, btype in oxy_neighbors:
            if btype.name == "DOUBLE":
                dbl_bond = True
            elif btype.name == "SINGLE":
                sgl_bond = True
        if not (dbl_bond and sgl_bond):
            continue

        # Now, count carbon neighbors (excluding the oxygens)
        carbon_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_neighbors.append(nbr)
        if len(carbon_neighbors) != 1:
            continue

        # We have a candidate terminal carboxyl carbon.
        candidate_found = True
        candidate_acid_atom = atom
        candidate_chain_neighbor = carbon_neighbors[0]
        break

    if not candidate_found:
        return False, "No terminal carboxyl group found (acid carbon not identified with required bonding)"

    # For the fatty acid chain to be valid, the acid carbon itself should not be in a ring.
    if candidate_acid_atom.IsInRing():
        return False, "Acid carbon is in a ring; not a linear fatty acid chain"

    # Now follow the chain from the unique carbon neighbor.
    # We require an unbranched, acyclic path. Here we do an iterative “walk”:
    chain_length = 1  # start with acid carbon
    current = candidate_chain_neighbor
    previous_idx = candidate_acid_atom.GetIdx()

    # Check: the chain neighbor should not itself be in a ring.
    if current.IsInRing():
        return False, "Chain carbon is in a ring; not a valid fatty acid chain"

    chain_length += 1  # counting the first chain carbon

    # Walk the chain until no unambiguous next carbon is found.
    while True:
        # For the current atom, get its carbon neighbors (excluding the one we came from)
        nbr_carbon_indices = []
        for nbr in current.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() == previous_idx:
                continue
            nbr_carbon_indices.append(nbr.GetIdx())
        
        # For a linear unbranched chain, there must be exactly one next carbon.
        if len(nbr_carbon_indices) == 0:
            # reached terminal end
            break
        if len(nbr_carbon_indices) > 1:
            return False, f"Chain is branched at carbon index {current.GetIdx()}; not a linear fatty acid chain"

        # Otherwise, continue along the unique path
        next_idx = nbr_carbon_indices[0]
        next_atom = mol.GetAtomWithIdx(next_idx)
        if next_atom.IsInRing():
            return False, f"Chain carbon at index {next_idx} is in a ring; not a valid fatty acid chain"
        chain_length += 1
        previous_idx = current.GetIdx()
        current = next_atom

    # Finally, check the chain length criterion (including the acid carbon).
    if chain_length < 13:
        return False, f"Chain length too short: {chain_length} carbons (< 13 required)"
    if chain_length > 22:
        return False, f"Chain length too long: {chain_length} carbons (> 22 allowed)"

    return True, f"Terminal carboxyl group found with linear chain of {chain_length} carbons"

# Example test cases (for local testing, uncomment if needed)
# test_smiles = [
#     ("O(O)[C@H](CCCCC)\\C=C\\CCCCCCCCCC(O)=O", "13R-HpOME(11E)"),
#     ("O[C@H]([C@@H](O)C/C=C/CCCCC)C(=O)CCCCCCC(O)=O", "(9R,10S,12Z)-9,10-Dihydroxy-8-oxo-12-octadecenoic acid"),
#     ("CCCCCCCCCCCCC(O)=O", "Tridecanoic acid"),
# ]
# for s, name in test_smiles:
#     result, reason = is_long_chain_fatty_acid(s)
#     print(f"SMILES: {s}\nNAME: {name}\nResult: {result}, Reason: {reason}\n")
    
# Note: This heuristic focuses on identifying a terminal carboxyl carbon that is bonded as expected
# and then “walking” the unbranched, acyclic carbon chain attached to it. Molecules with additional
# functional groups or cyclic/branched structures (such as carnitines or prostaglandins) will thus be 
# excluded from the classification. This approach may still require further refinement for edge cases.