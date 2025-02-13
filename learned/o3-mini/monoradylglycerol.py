"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: Monoradylglycerol
Definition:
    Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.
    
Our approach:
 1. Parse the SMILES.
 2. Look for a three-carbon linear fragment (all sp3 carbons connected by a single bond)
    that could serve as a glycerol backbone.
 3. For each candidate backbone, gather the oxygen substituents (neighbors that are not part of the chain).
    In “free” glycerol, each carbon has an –OH group; in a monoradylglycerol exactly one of those three oxygens will be substituted 
    (i.e. will have no attached hydrogen).
 4. Verify that exactly two oxygens are free (i.e. get at least one hydrogen) and exactly one oxygen is substituted (get zero hydrogens).
 5. For the substituted oxygen, follow its bond away from the backbone and (if attached to carbon) assess that there is at least 
    a modest carbon chain (length ≥ 6) present.
 6. If any candidate backbone passes these criteria, we return True with a reason.

If no candidate is found or the SMILES is invalid, we return False with an explanatory reason.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines whether a molecule is a monoradylglycerol.
    A monoradylglycerol is any glycerol that bears exactly one fatty (acyl, alkyl or alk-1-enyl) substituent
    (i.e. one of the three –OH groups is replaced) while the remaining two –OH groups are free.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as monoradylglycerol, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule.
    Chem.SanitizeMol(mol)
    
    # Helper function: recursively compute longest carbon chain length starting from a carbon atom.
    def longest_chain_length(atom, visited):
        # Only consider carbon atoms (atomic number 6).
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1  # current carbon counts as 1
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            # Only traverse through carbon atoms not yet visited.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                length = 1 + longest_chain_length(nbr, visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    # Get list of all carbon atoms that are sp3.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
    
    # We will look for any three-carbon chain (neighbors connected linearly) as candidate glycerol backbones.
    visited_chains = set()
    for a in carbons:
        for b in a.GetNeighbors():
            if b.GetAtomicNum() != 6 or b.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # a and b are connected
            for c in b.GetNeighbors():
                if c.GetAtomicNum() != 6 or c.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                # Ensure c is not a (avoid trivial 2-atom loops)
                if c.GetIdx() == a.GetIdx():
                    continue
                # We now have a chain: a - b - c. To avoid double counting, sort indices.
                chain = tuple(sorted([a.GetIdx(), b.GetIdx(), c.GetIdx()]))
                if chain in visited_chains:
                    continue
                visited_chains.add(chain)
                
                # Check if these three atoms form a connected path.
                # We require that a and c are terminal (each must be connected to b) and b is in the middle.
                if not (a in b.GetNeighbors() and c in b.GetNeighbors()):
                    continue
                
                # For each carbon in the candidate backbone, collect oxygen substituents (neighbors not in the backbone).
                backbone_idxs = {a.GetIdx(), b.GetIdx(), c.GetIdx()}
                oxygen_substituents = []
                for carbon in (a, b, c):
                    # For glycerol, we expect one oxygen attached to each carbon.
                    # We allow only oxygens (atomic number 8) that are not part of the backbone.
                    oxy_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_idxs]
                    if len(oxy_neighbors) != 1:
                        # We expect exactly one oxygen substituent per backbone carbon.
                        oxy_neighbors = []  # Invalidate candidate
                        break
                    oxygen_substituents.append((carbon, oxy_neighbors[0]))
                if len(oxygen_substituents) != 3:
                    continue
                
                # Count free hydroxyls vs substituted oxygen.
                free_count = 0
                substituted_idx = None
                for parent, oxy in oxygen_substituents:
                    # Using GetTotalNumHs gives the total count of (implicit+explicit) hydrogens.
                    # A free hydroxyl should have >= 1 hydrogen.
                    if oxy.GetTotalNumHs() > 0:
                        free_count += 1
                    else:
                        substituted_idx = oxy.GetIdx()
                # For a monoradylglycerol we need exactly 2 free hydroxyls and 1 substituted oxygen.
                if free_count != 2:
                    continue

                # Identify the substituted oxygen and check if it is linked to a fatty chain.
                substituted_atom = None
                for parent, oxy in oxygen_substituents:
                    if oxy.GetTotalNumHs() == 0:
                        substituted_atom = oxy
                        break
                if substituted_atom is None:
                    continue
                # Get neighbors of the substituted oxygen except the glycerol carbon.
                fatty_candidates = [nbr for nbr in substituted_atom.GetNeighbors() if nbr.GetIdx() not in backbone_idxs]
                if len(fatty_candidates) == 0:
                    continue
                # For each candidate that is carbon (or begins with a carbonyl carbon) check chain length.
                fatty_found = False
                for nbr in fatty_candidates:
                    # If the substituent oxygen is from an ester, the neighbor might be a carbonyl carbon.
                    # In that case, to evaluate chain length, look at the neighbor atoms that are carbons (excluding the oxygen).
                    if nbr.GetAtomicNum() == 6:
                        # For an ester, the pattern is O-C(=O)R. We check if nbr is carbonyl:
                        is_carbonyl = any(ne.GetAtomicNum() == 8 and (bond.GetBondTypeAsDouble() == 2) \
                                           for bond, ne in [(mol.GetBondBetweenAtoms(nbr.GetIdx(), x.GetIdx()), x) for x in nbr.GetNeighbors() if x.GetIdx() != substituted_atom.GetIdx()])
                        if is_carbonyl:
                            # For acyl chains, try to follow from the carbonyl carbon: find neighbor (other than the oxygen and glycerol)
                            for acyl_nbr in nbr.GetNeighbors():
                                if acyl_nbr.GetIdx() in (substituted_atom.GetIdx(),):
                                    continue
                                if acyl_nbr.GetAtomicNum() == 6:
                                    chain_len = longest_chain_length(acyl_nbr, set())
                                    if chain_len >= 6:
                                        fatty_found = True
                                        break
                        else:
                            # If not carbonyl then assume direct alkyl/alkenyl attachment.
                            chain_len = longest_chain_length(nbr, set())
                            if chain_len >= 6:
                                fatty_found = True
                    if fatty_found:
                        break
                if not fatty_found:
                    # The substituent did not appear to be a sufficiently long fatty chain.
                    continue

                # If we reach here, we have found a candidate glycerol backbone with exactly one substitution that appears to be fatty.
                return True, "Found glycerol backbone with exactly one fatty substituent (monoradylglycerol)"
    return False, "Could not identify a glycerol backbone with a single fatty substituent"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "O(C(=O)CCCCCCCCC/C=C/CCCCCC)CC(O)CO",  # MG(18:1(11E)/0:0/0:0)
        "CCCCCCCCCCCCCCCCCOC[C@@H](O)CO",         # 1-O-[(E)-hexadecen-1-yl]-sn-glycerol (simplified)
        "O(C(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO",  # 2-decanoylglycerol
        "CCCC(=O)OCC(O)CO"                        # monobutyrin (butyl is only 4 C so should fail the fatty chain check)
    ]
    for s in test_smiles:
        result, reason = is_monoradylglycerol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")