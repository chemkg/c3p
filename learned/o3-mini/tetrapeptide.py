"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino–acid residues connected by peptide linkages.
For a linear tetrapeptide the main backbone should contain exactly three peptide bonds in a continuous chain.
This improved version uses a simpler SMARTS pattern "N-C(=O)-C" (ignoring chirality) to capture a backbone peptide bond.
Then it examines connectivity by requiring that the alpha–carbon (the third atom from the match) of one bond is directly bonded
to the amide nitrogen (the first atom of the match) of the next, thereby identifying a continuous chain of three peptide bonds.
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if the given molecule (as a SMILES string) is a tetrapeptide.
    We use the strategy:
      1. Find all potential backbone peptide bonds via the SMARTS pattern "N-C(=O)-C".
         (This pattern ignores chirality and so is less strict but more robust to stereochemical specifications.)
      2. Build connectivity between these candidate peptide bonds: if the "C" (assumed alpha–carbon)
         of one candidate is directly bonded to the "N" of another candidate, they may be consecutively connected.
      3. Use depth–first search (DFS) to find a contiguous chain of three peptide bonds (thus 4 residues).
    
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
    # Using "N-C(=O)-C" rather than "N-C(=O)-[C]" to avoid issues with explicit chiral specifications.
    peptide_pattern = Chem.MolFromSmarts("N-C(=O)-C")
    if peptide_pattern is None:
        return False, "Could not construct peptide bond SMARTS pattern."
    
    # Get all substructure matches.
    # Each match returns a tuple of atom indices corresponding to (N, C, alpha–C)
    matches = mol.GetSubstructMatches(peptide_pattern)
    if not matches:
        return False, "No peptide bonds found."
    
    # Each match is considered a candidate backbone peptide bond.
    bonds = []
    for m in matches:
        bonds.append({'N': m[0], 'C': m[1], 'Ca': m[2]})
    
    # Build an adjacency mapping between candidate bonds.
    # We define that bond i -> bond j if the alpha–carbon (Ca) of bond i is directly bonded to the amide nitrogen (N) of bond j.
    adjacency = {i: [] for i in range(len(bonds))}
    for i in range(len(bonds)):
        for j in range(len(bonds)):
            if i == j:
                continue
            ca_i = bonds[i]['Ca']
            n_j = bonds[j]['N']
            if mol.GetBondBetweenAtoms(ca_i, n_j) is not None:
                adjacency[i].append(j)
    
    # Depth-first search (DFS) to find a continuous chain of 3 peptide bonds.
    # A chain of 3 peptide bonds corresponds to 4 residues.
    def dfs(current, depth, visited):
        # If we have chained 3 bonds then we have 4 connected residues.
        if depth == 3:
            return True
        for neighbor in adjacency[current]:
            if neighbor in visited:
                continue
            if dfs(neighbor, depth + 1, visited | {neighbor}):
                return True
        return False

    found_chain = False
    for i in range(len(bonds)):
        if dfs(i, 1, {i}):
            found_chain = True
            break
    
    if found_chain:
        return True, "Molecule contains a continuous backbone chain with 3 peptide bonds (linking 4 residues) consistent with a tetrapeptide."
    else:
        return False, "Did not find a continuous chain of 3 backbone peptide bonds connecting 4 residues; peptide bond connectivity does not match a tetrapeptide."

# Example usage (can be removed or commented out in production):
if __name__ == "__main__":
    # Glu-Lys-Trp-Ala example tetrapeptide:
    sample_smiles = "C[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](N)CCC(O)=O)C(O)=O"
    result, reason = is_tetrapeptide(sample_smiles)
    print(result, reason)