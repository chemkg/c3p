"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
#!/usr/bin/env python3
"""
Classifies: fatty acid methyl ester (FAME)

A FAME is defined as a fatty acid ester obtained from the formal 
condensation of a fatty acid with methanol. Our approach is:
  1. Look for an ester group matching the SMARTS [C:1](=[O:2])[O:3][C:4].
  2. Check that the alkoxy fragment ([C:4]) is a methyl group.
  3. Ensure that the carbonyl carbon (atom :1) carries exactly one carbon substituent 
     (the acyl chain) aside from the ester oxygen.
  4. “Extract” the acyl fragment (by a DFS that forbids crossing the bond to the ester oxygen)
     and require that:
       • the acyl fragment contains at least 4 carbon atoms (including the carbonyl carbon),
       • the acyl fragment does not contain any ring atoms,
       • and it does not contain an additional carboxyl group.
These added criteria help to filter out esters that are merely methoxy‐substituted 
or part of more complex scaffolds.
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester (FAME) based on its SMILES.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if defined as a FAME, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define the ester SMARTS: 
    # [C:1](=[O:2])[O:3][C:4] 
    ester_smarts = "[C:1](=[O:2])[O:3][C:4]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "No ester functional group matching [C](=[O])[O]C found"
    
    # Helper: check if an atom is a methyl group (only one heavy neighbor)
    def is_methyl(atom):
        if atom.GetAtomicNum() != 6:
            return False
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        return len(heavy_neighbors) == 1

    # Helper: perform a DFS from a given start atom, but do not cross
    # the bond between carbonyl and the ester oxygen (to avoid the methoxy group).
    def dfs_acyl(atom, forbidden_bond, visited):
        visited.add(atom.GetIdx())
        atoms = {atom.GetIdx()}
        for bond in atom.GetBonds():
            # skip the forbidden bond
            if bond.GetIdx() == forbidden_bond:
                continue
            # get the other atom on this bond
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetIdx() in visited:
                continue
            # continue the DFS unconditionally (we want the whole connected fragment)
            atoms.update(dfs_acyl(nbr, forbidden_bond, visited))
        return atoms

    # Define a SMARTS to detect a carboxyl group pattern within the acyl fragment.
    # This simple pattern [C](=O)[O] should match a carboxyl group.
    carboxyl_smarts = Chem.MolFromSmarts("C(=O)[O]")
    
    # Loop over each ester match candidate
    for match in matches:
        # In the match: 
        # match[0] -> carbonyl carbon, match[1] -> carbonyl oxygen (double-bonded),
        # match[2] -> ester oxygen, match[3] -> alkoxy carbon.
        carbonyl = mol.GetAtomWithIdx(match[0])
        ester_oxygen = mol.GetAtomWithIdx(match[2])
        alkoxy = mol.GetAtomWithIdx(match[3])
        
        # Check that the alkoxy group is a methyl group.
        if not is_methyl(alkoxy):
            continue  # Not methoxy; try next match
        
        # Identify the acyl chain:
        # It should be the unique carbon neighbor of the carbonyl (other than the ester oxygen)
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors() 
                          if nbr.GetIdx() != ester_oxygen.GetIdx() and nbr.GetAtomicNum() == 6]
        if len(acyl_neighbors) != 1:
            continue  # ambiguous or missing acyl chain
        acyl_start = acyl_neighbors[0]
        
        # (Optionally) if carbonyl or acyl_start are in a ring, reject this match.
        if carbonyl.IsInRing() or acyl_start.IsInRing():
            continue
        
        # Identify the bond between carbonyl and the ester oxygen (to prevent crossing into the methoxy part)
        bond_to_block = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), ester_oxygen.GetIdx())
        if bond_to_block is None:
            continue
        forbidden_bond_idx = bond_to_block.GetIdx()
        
        # Use DFS starting from the carbonyl to get the entire acyl fragment,
        # but do not cross the bond to the ester oxygen.
        acyl_atom_indices = dfs_acyl(carbonyl, forbidden_bond_idx, set())
        
        # Count number of carbon atoms in the acyl fragment.
        n_carbons = sum(1 for idx in acyl_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbons < 4:
            continue  # Acyl portion too short to be a typical fatty acid
        
        # Check that none of the atoms in the acyl fragment is in a ring.
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in acyl_atom_indices):
            continue  # fatty acids are expected to have an acyclic acyl chain
        
        # Build a submol for the acyl fragment (for additional analysis).
        # Note: Chem.PathToSubmol accepts a list of atom indices.
        submol = Chem.PathToSubmol(mol, list(acyl_atom_indices))
        # Count the number of carboxyl groups (pattern "C(=O)[O]") in the acyl fragment.
        cps = submol.GetSubstructMatches(carboxyl_smarts)
        if len(cps) > 1:
            # More than one carboxyl group suggests a diacid ester rather than a FAME
            continue
        
        # If we reach here, we have found an ester match whose alkoxy part is methyl,
        # the carbonyl carbon bears a unique acyl chain that is acyclic, sufficiently long,
        # and shows only a single carboxyl functionality.
        return True, "Contains a methoxy ester group with a proper acyl chain (FAME identified)"
    
    return False, "Ester functional group found but does not match fatty acid methyl ester criteria"


# (Optional) Quick testing when module is run directly.
if __name__ == "__main__":
    # Test with an example that should be classified as FAME (methyl octanoate)
    test_smiles = "O=C(OC)CCCCCCCC"
    result, reason = is_fatty_acid_methyl_ester(test_smiles)
    print(result, reason)