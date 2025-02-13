"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: Tetrapeptide
Definition: A tetrapeptide is any molecule that contains four amino‐acid residues connected by peptide linkages.
For a linear tetrapeptide the main backbone should contain exactly three peptide bonds in a continuous chain.
This program uses a two–step approach:
  1. It searches for backbone “peptide bond” fragments using a SMARTS pattern.
     (The pattern is designed to capture an amide bond along with the next (alpha) carbon.)
  2. It then looks for a contiguous chain of three peptide bonds by requiring that the “alpha-carbon”
     (from one match) is directly bonded to the amide nitrogen of the next match.
This connectivity criterion should weed out extra (side–chain or protecting group) amide bonds.
"""

from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if the given molecule (as a SMILES string) is a tetrapeptide.
    The strategy is:
      1. Use a SMARTS pattern that roughly captures a backbone amide bond together with its subsequent alpha–carbon.
         In a linear peptide the backbone (ignoring terminal modifications) is, in simplified form:
           H2N–CH(R)–C(=O)–NH–CH(R)–C(=O)–...–C(=O)–[rest]
         Here we use a pattern "N–C(=O)–[C]" where:
           • the N is (at least partially) in an amide bond,
           • the C(=O) is the carbonyl carbon,
           • and the [C] is assumed to be the alpha–carbon of the residue on the C–terminal side.
      2. From all matches of the peptide bond pattern, we then try to “chain” them.
         In a proper tetrapeptide, there is one contiguous chain of exactly three backbone bonds connecting four residues.
         We say that one peptide bond match (with atoms (N, C, Ca)) is followed by another if the alpha–carbon (Ca)
         of the first bond is directly connected (by a bond) to the amide nitrogen (N) of the second.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a tetrapeptide, False otherwise.
      str: Explanation for the result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to capture a potential backbone peptide bond.
    # The SMARTS "N-C(=O)-[C]" means:
    #   an N atom directly bonded to a carbon that is double-bonded to O,
    #   and that same carbon is bonded to a carbon (assumed to be the alpha–carbon).
    #
    # Note: This pattern may also pick up extra amide bonds (from side–chains or modifications).
    peptide_pattern = Chem.MolFromSmarts("N-C(=O)-[C]")
    if peptide_pattern is None:
        return False, "Could not construct peptide bond SMARTS pattern."
    
    # Get all substructure matches.
    # Each match is a tuple of atom indices (n_idx, c_idx, ca_idx)
    matches = mol.GetSubstructMatches(peptide_pattern)
    if not matches:
        return False, "No peptide bonds found."
    
    # Each match is taken to represent a backbone peptide bond.
    # However, extra amide bonds (not along the main chain) may also be found.
    # We now try to look for a continuous chain of peptide bonds.
    # Our idea: if one match (bond1) with atoms (N1, C1, Ca1)
    # and another match (bond2) with atoms (N2, C2, Ca2) are “consecutive”
    # along the main chain then the alpha–carbon Ca1 should be bonded to the amide nitrogen N2.
    # We store each match as a dictionary for easier access.
    bonds = []
    for m in matches:
        bonds.append({'N': m[0], 'C': m[1], 'Ca': m[2]})
    
    # Build a mapping (graph) between matches.
    # For two matches bond_i and bond_j, we say i --> j if the alpha–carbon of bond_i is directly bonded
    # to the amide nitrogen of bond_j.
    # This check uses the molecule’s bond connectivity.
    adjacency = {i: [] for i in range(len(bonds))}
    for i in range(len(bonds)):
        for j in range(len(bonds)):
            if i == j:
                continue
            # Check if the alpha-carbon of bond i is directly bonded to the nitrogen of bond j.
            ca_i = bonds[i]['Ca']
            n_j = bonds[j]['N']
            if mol.GetBondBetweenAtoms(ca_i, n_j) is not None:
                adjacency[i].append(j)
    
    # Now we search for a chain (path) of consecutive bonds of length 3.
    # A path length of 3 in terms of peptide bonds means 4 residues connected by 3 backbone peptide bonds.
    # We perform a DFS from each bond.
    def dfs(start, depth, visited):
        # depth counts how many bonds are in the chain so far.
        if depth == 3:
            return True  # found 3 consecutive bonds
        for neighbor in adjacency[start]:
            # To avoid cycles we do not re-visit the same match in this chain.
            if neighbor in visited:
                continue
            # Extend the chain
            if dfs(neighbor, depth + 1, visited | {neighbor}):
                return True
        return False
    
    found_chain = False
    for i in range(len(bonds)):
        if dfs(i, 1, {i}):
            found_chain = True
            break
    
    if found_chain:
        return True, "Molecule contains a continuous backbone chain with 3 peptide bonds (linking 4 residues) consistent with a tetrapeptide"
    else:
        return False, "Did not find a continuous chain of 3 backbone peptide bonds connecting 4 residues; peptide bond pattern count or connectivity does not match a tetrapeptide"

# Example usage:
if __name__ == "__main__":
    # An example tetrapeptide (Glu-Lys-Trp-Ala)
    sample_smiles = "C[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](N)CCC(O)=O)C(O)=O"
    result, reason = is_tetrapeptide(sample_smiles)
    print(result, reason)