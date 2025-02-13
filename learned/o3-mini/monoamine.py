"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine
Definition: An aralylamino compound which contains one amino group connected 
to an aromatic ring by a two‐carbon chain.
"""

from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an aralylamino compound that contains one amino 
    group connected to an aromatic ring by a two‐carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as a monoamine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -- Exclude molecules that are likely amino acid derivatives --
    # Typical amino acids contain a carboxyl group. Here, we look for a pattern that has:
    # a carbonyl (C=O) and a hydroxyl (O) attached to that carbon.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if carboxyl_pattern and mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Contains carboxyl group, likely an amino acid derivative rather than a monoamine"
    
    # -- Define a refined SMARTS pattern for the aromatic ethylamine chain --
    # The pattern looks for:
    #  1. An aromatic atom [c] (which we will later check to be in a six‐membered ring),
    #  2. Followed by two aliphatic CH2 groups, and then 
    #  3. A nitrogen (which can be primary, secondary, or even protonated).
    # We label the atoms so we can later check the aromatic atom and the nitrogen.
    smarts = "[c:1][CH2:2][CH2:3][N:4]"
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        return False, "Error in SMARTS pattern definition"
    
    # Find all substructure matches for our pattern.
    matches = mol.GetSubstructMatches(pattern, uniquify=True)
    
    # Filter matches:
    # We require that the aromatic atom (match[0]) is in a ring of size 6
    # and that the two chain carbons (match[1] and match[2]) are not aromatic.
    valid_matches = []
    rings = mol.GetRingInfo().AtomRings()
    for match in matches:
        aromatic_atom = mol.GetAtomWithIdx(match[0])
        # Check if this atom is in a ring of exactly 6 atoms (a benzene-type ring)
        in_6_membered_ring = any(len(ring) == 6 and match[0] in ring for ring in rings)
        if not in_6_membered_ring:
            continue
        # Ensure the two carbon atoms in the chain are aliphatic.
        atom2 = mol.GetAtomWithIdx(match[1])
        atom3 = mol.GetAtomWithIdx(match[2])
        if atom2.GetIsAromatic() or atom3.GetIsAromatic():
            continue
        valid_matches.append(match)
    
    if not valid_matches:
        return False, "No valid aromatic ethylamine chain found"

    # Count the number of unique nitrogen atoms identified in the valid matches.
    amine_nitrogens = set(match[3] for match in valid_matches)
    
    if len(amine_nitrogens) != 1:
        return False, f"Found {len(amine_nitrogens)} distinct aromatic ethylamine chain(s); expected exactly one"
    
    # Passed all tests: exactly one valid aryl–ethylamine chain was found and no disqualifying groups were detected.
    return True, "Found a single aromatic ethylamine chain (aralylamino group) indicative of a monoamine"

# Example usage (to test, uncomment and run):
# smiles_examples = [
#     "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",  # (S)-dobutamine
#     "C1=C(C(=CC(=C1O)O)[N+](=O)[O-])CCN",          # 4-(2-aminoethyl)-5-nitrobenzene-1,2-diol
#     "NCCc1ccc(O)c(O)c1",                           # dopamine hydrochloride
#     "C[C@H](N)[C@H](O)c1ccc(O)c(O)c1"              # (-)-alpha-Methylnoradrenaline
# ]
# for smi in smiles_examples:
#     result, reason = is_monoamine(smi)
#     print(f"SMILES: {smi}\nMonoamine: {result}\nReason: {reason}\n")