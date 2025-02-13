"""
Classifies: CHEBI:61384 sulfolipid
"""
#!/usr/bin/env python
"""
Classifies: sulfolipid
Definition: A compound containing a sulfonic acid residue (–SO3H or –SO3–) that is joined (directly or via a bridging oxygen) by a carbon–sulfur or carbon–oxygen–sulfur linkage to a lipid.
Heuristic strategy:
  1. Parse the SMILES.
  2. To avoid false positives on simple sulfonic acids (e.g. hexadecane-1-sulfonic acid), we require that the overall molecule shows “complexity”
     (for example, at least one ring or one amide bond).
  3. For every sulfur atom, check if it is part of a sulfonate/sulfonic acid subgroup by requiring at least three oxygen neighbors (with ≥2 double bonds).
  4. For each such S, check if one of its neighbor atoms (directly or via an oxygen bridge) is a carbon that is “aliphatic” (sp3, non‐aromatic)
     and is part of a contiguous chain of at least 10 sp3 carbons.
     The chain search (via DFS) is restricted to non‐aromatic sp3 carbons connected by single bonds.
  5. If a candidate is found, return True together with an explanation.
  
If no candidate S is found that is connected to a sufficiently long lipid chain, return False.
  
Note: This algorithm is heuristic and may mis‐categorize borderline structures.
"""

from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule should be classified as a sulfolipid based on a heuristic:
      (a) the molecule must be “complex” (contain a ring or an amide bond), and
      (b) it must contain at least one sulfonate/sulfonic acid-like moiety (S bound to ≥3 O atoms, with ≥2 double bonds), and
      (c) one neighbor (directly or via an oxygen) of that S is a carbon that is part of a contiguous aliphatic chain of at least 10 sp3 carbons.
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): A tuple where the first element is True if the molecule is classified as a sulfolipid,
                   False otherwise. The second element is a textual explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # First check for minimal molecular complexity (e.g. rings or amide bonds)
    if not mol.GetRingInfo().NumRings() and not mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)N")):
        # If the molecule is acyclic and does not contain an amide bond, then it might be a simple sulfonic acid.
        return False, "Molecule appears too simple (acyclic and lacks amide bonds) to be a sulfolipid."
    
    # DFS helper to walk a contiguous chain of aliphatic sp3 carbons (single bonds only)
    def dfs_chain(atom, visited):
        # Count the current atom and search for neighbors that are carbon, sp3 and non‐aromatic.
        length = 1
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Only traverse if neighbor is carbon (atomic number 6), non‐aromatic and the bond is single.
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and bond.GetBondType() == Chem.BondType.SINGLE:
                if nbr.GetIdx() in visited:
                    continue
                new_visited = visited.copy()
                new_visited.add(nbr.GetIdx())
                branch_length = 1 + dfs_chain(nbr, new_visited)
                if branch_length > length:
                    length = branch_length
        return length

    # Now loop over all sulfur atoms as potential sulfonate/sulfonic acid centers.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # not S
        
        # For candidate S, gather oxygen neighbors and count bond orders.
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                oxy_neighbors.append(nbr)
        # A sulfonate/sulfonic acid ideally has at least three O substituents
        if len(oxy_neighbors) < 3:
            continue
        
        double_bond_count = 0
        single_bond_count = 0
        for nbr in oxy_neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                double_bond_count += 1
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                single_bond_count += 1
        if double_bond_count < 2 or single_bond_count < 1:
            continue   # does not fulfill a typical sulfonic acid pattern

        checked_candidates = False
        # Now check all neighbors of S.
        # We allow the possibility that S is connected directly to a carbon with a lipid chain...
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                chain_len = dfs_chain(nbr, {nbr.GetIdx()})
                if chain_len >= 10:
                    explanation = (f"Contains a sulfonate group (S(=O)(=O)(O) or S(=O)(=O)[O-]) directly attached via a C–S bond "
                                   f"to an aliphatic chain of {chain_len} contiguous sp3 carbons.")
                    return True, explanation
                else:
                    checked_candidates = True
            # Alternatively, the S may be attached via an intervening oxygen.
            if nbr.GetAtomicNum() == 8:
                # For each oxygen neighbor, check if it itself is bound to a carbon (other than S) that yields a long lipid chain.
                for o_nbr in nbr.GetNeighbors():
                    if o_nbr.GetIdx() == atom.GetIdx():
                        continue
                    if o_nbr.GetAtomicNum() == 6 and not o_nbr.GetIsAromatic():
                        chain_len = dfs_chain(o_nbr, {o_nbr.GetIdx()})
                        if chain_len >= 10:
                            explanation = (f"Contains a sulfonate group (S(=O)(=O)(O) or S(=O)(=O)[O-]) attached via an O bridge "
                                           f"to a carbon that is part of an aliphatic chain of {chain_len} contiguous sp3 carbons.")
                            return True, explanation
                        else:
                            checked_candidates = True

        if checked_candidates:
            # We found candidate neighbor(s) but their lipid chain length was too short.
            return False, "Candidate sulfonate observed but the attached aliphatic chain is too short (< 10 contiguous sp3 carbons)."
    
    return False, "No sulfonate (sulfonic acid) group attached (directly or via an oxygen) to a sufficiently long aliphatic chain was found."


# For testing, you might uncomment one or more examples below:
# test_smiles = "[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(CO[C@@H]1O[C@@H]([C@@H]([C@@H]([C@H]1O)OS(=O)(=O)O)O)CO)NC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O"
# result, reason = is_sulfolipid(test_smiles)
# print(result, reason)

# If uncertain, one is allowed to return (None, None) in borderline cases.