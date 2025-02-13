"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N-acylglycine
Definition: An N-acyl-amino acid in which the amino acid specified is glycine.
The key motif is: R-C(=O)-N-CH2-C(=O)O.
This implementation uses a relaxed SMARTS pattern and then post-filters
matches to ensure that the glycine α-carbon has exactly two hydrogen atoms
and that the amide nitrogen is terminal (only two heavy-atom bonds).
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine must contain an acyl group attached to glycine via an amide bond,
    resulting in the motif: R-C(=O)-N-CH2-C(=O)O.
    To avoid false positives (e.g. when a glycine residue is embedded in a peptide),
    we check that the glycine unit is exactly represented (the α-carbon has two hydrogens)
    and that the amide nitrogen is connected only to the acyl carbon and the glycine α-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a valid N-acylglycine motif, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a relaxed SMARTS pattern for the motif.
    # We capture four atoms with explicit mapping:
    #   [C:1](=O)   --> acyl carbonyl (part of the R-C(=O) group)
    #   [N:2]       --> amide nitrogen
    #   [C:3]       --> glycine α-carbon (will be post-filtered for exactly 2 H's)
    #   [C:4](=O)[O] --> carboxylic acid carbon in glycine ([O] can be -OH or O-)
    #
    # This pattern is intentionally permissive.
    motif_smarts = "[C:1](=O)-[N:2]-[C:3]-[C:4](=O)[O]"
    pattern = Chem.MolFromSmarts(motif_smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern."
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "The molecule does not contain the N-acylglycine moiety (R-C(=O)-N-CH2-C(=O)O)."
    
    # Iterate over each match candidate and further enforce that:
    #   - The glycine α-carbon (mapped as atom 3) has exactly two hydrogen atoms.
    #   - The amide nitrogen (mapped as atom 2) has exactly two heavy-atom neighbors.
    for match in matches:
        # Retrieve atoms by their mapped indices.
        acylC = mol.GetAtomWithIdx(match[0])
        amideN = mol.GetAtomWithIdx(match[1])
        alphaC = mol.GetAtomWithIdx(match[2])
        acidC = mol.GetAtomWithIdx(match[3])
        
        # Check that the alpha (glycine) carbon has exactly two hydrogens.
        # This distinguishes glycine (CH2) from other amino acids.
        # Note: GetTotalNumHs() returns both implicit and explicit hydrogens.
        if alphaC.GetTotalNumHs() != 2:
            continue  # not glycine
        
        # Check that the amide nitrogen (the one connecting acyl and glycine)
        # is only connected to two heavy atoms (the acyl carbon and the glycine carbon).
        # Count neighbors that are not hydrogen:
        n_heavy_neighbors = sum(1 for nbr in amideN.GetNeighbors() if nbr.GetAtomicNum() > 1)
        if n_heavy_neighbors != 2:
            continue  # the N is further substituted (likely part of a peptide chain)
        
        # If this match passes the filters, we assume a valid N-acylglycine motif was found.
        return True, "Molecule contains the N-acylglycine moiety (R-C(=O)-N-CH2-C(=O)O)."
    
    # If no match passes the additional checks, return False.
    return False, "The molecule does not contain a valid N-acylglycine moiety after post-filtering."

# An example test can be run by executing this file directly.
if __name__ == "__main__":
    # Use N-benzoylglycine as an example.
    test_smiles = "OC(=O)CNC(=O)c1ccccc1"
    result, reason = is_N_acylglycine(test_smiles)
    print(f"SMILES: {test_smiles}\nClassification: {result}\nReason: {reason}\n")