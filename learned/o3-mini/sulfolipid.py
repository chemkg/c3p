"""
Classifies: CHEBI:61384 sulfolipid
"""
#!/usr/bin/env python
"""
Classifies: sulfolipid
Definition: A compound containing a sulfonic acid residue joined by a carbon–sulfur bond to a lipid 
            (where “lipid” is defined as a contiguous chain of at least 10 aliphatic [sp3, non‐aromatic] carbons).

Strategy:
  1. Parse the SMILES.
  2. Identify sulfur atoms that are part of a sulfonic acid group.
     We require that the sulfur is bonded to at least three oxygens – two via double bonds and one via a single bond.
  3. For each such S, check if at least one neighbor is a carbon (indicating a C–S bond).
  4. For each candidate carbon, use DFS to find the maximum length of contiguous aliphatic (sp3 non‐aromatic) carbon chain (via single bonds).
  5. Return True if any chain is at least 10 carbons long.
  
If either no sulfonic acid group or no such lipid chain is found, return False.
"""

from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid by verifying that:
      (a) it contains a sulfonic acid (or sulfonate) residue (S(=O)(=O)(O) or S(=O)(=O)[O-]),
      (b) that sulfonate sulfur is directly bonded to a carbon (i.e., via a C–S bond),
      (c) and that the carbon is part of a contiguous aliphatic chain of at least 10 sp3 carbons.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a sulfolipid, otherwise False.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Identify candidate sulfonic acid/sulfonate groups:
    # For each sulfur atom in the molecule, we look at its neighbors.
    candidate_found = False
    candidate_explanation = ""
    
    # Helper DFS to compute the longest contiguous chain of sp3, non-aromatic carbons (starting from an atom)
    def dfs(atom, visited):
        max_length = 1  # current atom is counted
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # walk only through carbon atoms that are sp3 (non-aromatic) and connected by a SINGLE bond
            if (nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic()
                and bond.GetBondType() == Chem.BondType.SINGLE
                and nbr.GetIdx() not in visited):
                new_visited = visited.copy()
                new_visited.add(nbr.GetIdx())
                chain_length = 1 + dfs(nbr, new_visited)
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # Loop over all sulfur atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16: # not sulfur
            continue

        # Check neighbors of S for oxygens.
        oxy_neighbors = []
        other_neighbors = []  # atoms other than oxygen
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                oxy_neighbors.append(nbr)
            else:
                other_neighbors.append(nbr)
        
        # We expect a sulfonic acid group to have at least 3 oxygen neighbors.
        # And among the oxygens, ideally 2 are connected by double bonds.
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
        # Require at least two double bonds and one single bond.
        if double_bond_count < 2 or single_bond_count < 1:
            continue

        # Now check among the other neighbors (or even oxygen neighbors) if there is a carbon directly attached.
        # We are looking for a C–S bond.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                # We have a candidate C–S bond.
                # Now get the contiguous aliphatic carbon chain (starting from that carbon).
                # To be conservative, check that the candidate carbon is sp3 (non-aromatic)
                if nbr.GetIsAromatic():
                    continue
                chain_length = dfs(nbr, {nbr.GetIdx()})
                if chain_length >= 10:
                    candidate_explanation = (f"Contains sulfonic acid group (S(=O)(=O)(O) or S(=O)(=O)[O-]) "
                                             f"joined via a C–S bond to a lipid chain of {chain_length} contiguous sp3 carbons.")
                    return True, candidate_explanation
                else:
                    # For debugging purposes record the chain length for this candidate if too short.
                    candidate_explanation = (f"Candidate C–S bond found but the attached aliphatic chain is only {chain_length} carbons long.")
                    candidate_found = True
        # End of sulfur loop

    if candidate_found:
        return False, candidate_explanation
    else:
        return False, "No sulfonic acid group attached via a C–S bond to a sufficiently long aliphatic (lipid) chain was found."

# For testing (you can uncomment and test with example SMILES):
# test_smiles = "[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(CO[C@@H]1O[C@@H]([C@@H]([C@@H]([C@H]1O)OS(=O)(=O)O)O)CO)NC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O"
# result, reason = is_sulfolipid(test_smiles)
# print(result, reason)