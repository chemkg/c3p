"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine
Definition: An aralylamino compound that contains one amino group connected to an aromatic ring by a two‐carbon chain.
Examples of structures belonging to this class include dopamine, tyramine, adrenaline, and others.
"""

from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an aralylamino compound that contains one amino 
    group connected to an aromatic ring by a two-carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as a monoamine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an aromatic ethylamine chain:
    # "[c]CC[N]" means: an aromatic carbon ([c]) connected to two aliphatic carbons (C C)
    # followed by a nitrogen atom ([N]). This roughly captures the aryl-alkylamine motif.
    pattern = Chem.MolFromSmarts("[c]CC[N]")
    if pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Find all substructure matches for our pattern in the molecule.
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No aryl-ethylamine chain (aralylamino functionality) found"
    
    # Because of symmetry, the same amino nitrogen may appear in several matches.
    # Hence, we collect the unique indices of the nitrogen atoms (position 3 in the match tuple).
    amino_nitrogens = set(match[3] for match in matches)
    if len(amino_nitrogens) == 0:
        return False, "No amino group in the aryl-ethylamine pattern detected"
    if len(amino_nitrogens) > 1:
        return False, f"Multiple ({len(amino_nitrogens)}) aryl-ethylamine chains found; expected one"
    
    # Optionally, check that the connecting chain contains exactly two carbon atoms.
    # Here we verify that the atoms in positions 1 and 2 of the match are indeed aliphatic carbons.
    match = matches[0]
    atom1 = mol.GetAtomWithIdx(match[1])
    atom2 = mol.GetAtomWithIdx(match[2])
    if atom1.GetIsAromatic() or atom2.GetIsAromatic():
        return False, "The connecting chain does not consist solely of aliphatic carbons"
    
    return True, "Found a single aromatic ethylamine chain (aralylamino group) indicative of a monoamine"

# Example usage (uncomment below lines to test):
# smiles_examples = [
#     "NCCc1ccc(O)c(O)c1",       # dopamine
#     "NCCc1ccc(O)cc1",          # tyramine
#     "CNC[C@@H](O)c1ccc(O)c(O)c1" # (S)-adrenaline
# ]
# for smi in smiles_examples:
#     result, reason = is_monoamine(smi)
#     print(f"SMILES: {smi}\nClassified as monoamine: {result}\nReason: {reason}\n")