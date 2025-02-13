"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N-acylglycine
Definition: An N-acyl-amino acid in which the amino acid specified is glycine.
The key motif is: R-C(=O)-N-CH2-C(=O)O.
This implementation uses an initial SMARTS pattern that is relaxed and then post-filters
matches to enforce that the glycine α‐carbon has exactly two hydrogens and that
the amide nitrogen (bridging the acyl group and the glycine unit) is terminal.
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine must contain an acyl group attached to glycine via an amide bond,
    resulting in the motif: R-C(=O)-N-CH2-C(=O)O. To avoid false positives when glycine
    residues are embedded within peptides, we check that (a) the glycine alpha-carbon
    has exactly two hydrogens (i.e. it is CH2) and (b) that the amide nitrogen only bears
    connections to the acyl group and the glycine alpha-carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a valid N-acylglycine motif is detected, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens so that the hydrogen counts in the glycine CH2 unit are explicit.
    mol = Chem.AddHs(mol)
    
    # Define a relaxed SMARTS pattern for the key motif.
    # Here we define four atoms with mappings:
    #   [C:1](=O)   --> acyl carbon (from the R-C(=O) group)
    #   [N:2]       --> amide nitrogen linking acyl to glycine
    #   [C:3]       --> glycine alpha-carbon (which must be CH2)
    #   [C:4](=O)[O] --> carboxylic acid carbon (glycine C-terminus)
    # Note: the pattern is written with explicit single bonds ('-').
    motif_smarts = "[C:1](=O)-[N:2]-[C:3]-[C:4](=O)[O]"
    pattern = Chem.MolFromSmarts(motif_smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern."
    
    # Search for substructure matches using the pattern.
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "The molecule does not contain the N-acylglycine motif (R-C(=O)-N-CH2-C(=O)O)."
    
    # For each match candidate, check:
    #   (1) the glycine alpha-carbon (mapped as atom 3) has exactly 2 hydrogens.
    #   (2) the amide nitrogen (mapped as atom 2) is terminal (only two heavy-atom neighbors).
    for match in matches:
        # Map each component of the motif.
        acylC = mol.GetAtomWithIdx(match[0])   # acyl carbonyl carbon (from the R-C(=O) group)
        amideN = mol.GetAtomWithIdx(match[1])    # amide nitrogen bridging groups
        alphaC = mol.GetAtomWithIdx(match[2])    # glycine α‐carbon, expected to be CH2
        acidC = mol.GetAtomWithIdx(match[3])     # glycine carboxyl carbon; no further filtering needed here
        
        # Ensure that the glycine alpha carbon (atom mapped as 3) has exactly two hydrogens.
        if alphaC.GetTotalNumHs() != 2:
            continue  # Not glycine; try next match.
        
        # Check that the amide nitrogen is only bonded to two heavy atoms.
        # (It should be connected only to the acyl carbon and the alpha carbon.)
        heavy_neighbors = [nbr for nbr in amideN.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 2:
            continue  # The amide nitrogen has extra substituents.
        
        # If both checks pass, we consider the N-acylglycine motif to be present.
        return True, "Molecule contains the N-acylglycine motif (R-C(=O)-N-CH2-C(=O)O)."
    
    # If no candidate match passes the additional filtering, return False.
    return False, "The molecule does not contain a valid N-acylglycine moiety after filtering."

# Example test cases.
if __name__ == "__main__":
    # Example: N-benzoylglycine.
    test_smiles = "OC(=O)CNC(=O)c1ccccc1"
    result, reason = is_N_acylglycine(test_smiles)
    print(f"SMILES: {test_smiles}\nClassification: {result}\nReason: {reason}\n")
    
    # You can test with other examples (such as N-acetylglycine and N-cinnamoylglycine) below.
    examples = [
        ("N-acetylglycine", "CC(=O)NCC(O)=O"),
        ("N-cinnamoylglycine", "OC(=O)CNC(=O)\\C=C\\c1ccccc1")
    ]
    for name, smiles in examples:
        res, rsn = is_N_acylglycine(smiles)
        print(f"Name: {name}\nSMILES: {smiles}\nClassification: {res}\nReason: {rsn}\n")