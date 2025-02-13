"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine
Definition: An aralylamino compound which contains one amino group connected 
to an aromatic ring by a two‐carbon chain. Monoamines are derived from aromatic
amino acids via decarboxylation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an aralylamino compound that contains one amino 
    group connected to an aromatic ring by a two‐carbon chain. In catecholamines
    the first carbon may carry a hydroxyl substituent. Additionally, typical monoamines
    are relatively small (here we require a molecular weight of 500 Da or less) and 
    must not have a carboxyl group (to avoid amino acid derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a monoamine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a carboxyl group (typical for amino acid derivatives).
    # This SMARTS searches for a carbonyl (C=O) with an -OH or -O^- attached.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if carboxyl_pattern and mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Contains carboxyl group, likely an amino acid derivative rather than a monoamine"
    
    # Prepare two SMARTS patterns:
    # Pattern 1: A simple aromatic ethylamine chain: aromatic carbon – CH2 – CH2 – nitrogen.
    pattern1 = Chem.MolFromSmarts("[c:1][CH2:2][CH2:3][N:4]")
    # Pattern 2: A variant found in catecholamines where the first carbon carries an -OH:
    # aromatic carbon – CH(OH) – CH2 – nitrogen.
    pattern2 = Chem.MolFromSmarts("[c:1][C;!a;H1]([OX2H])[CH2:3][N:4]")
    
    # This list will hold valid chain matches (as tuples of atom indices).
    valid_matches = []
    # Get ring information for filtering aromatic atoms in 6-membered rings.
    rings = mol.GetRingInfo().AtomRings()
    
    for pattern in (pattern1, pattern2):
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern, uniquify=True)
        for match in matches:
            # match indices: 
            #   match[0]: aromatic atom from the ring
            #   match[1] and match[2]: the two linker carbons
            #   match[3]: the nitrogen
            aromatic_idx = match[0]
            # Check that the aromatic atom is in at least one six-membered ring
            in_6_membered_ring = any(len(ring) == 6 and aromatic_idx in ring for ring in rings)
            if not in_6_membered_ring:
                continue
            # Check that the two linker carbons are not aromatic.
            atom1 = mol.GetAtomWithIdx(match[1])
            atom2 = mol.GetAtomWithIdx(match[2])
            if atom1.GetIsAromatic() or atom2.GetIsAromatic():
                continue
            valid_matches.append(match)
    
    if not valid_matches:
        return False, "No valid aromatic ethylamine chain found"
    
    # Count unique nitrogen atoms from all valid matches
    amine_nitrogens = set(match[3] for match in valid_matches)
    if len(amine_nitrogens) != 1:
        return False, f"Found {len(amine_nitrogens)} distinct aromatic ethylamine chain(s); expected exactly one"

    # Use molecular weight as an additional filter; monoamines are usually small.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} Da is too high for a typical monoamine"
    
    # If all checks pass, the molecule is classified as a monoamine.
    return True, "Found a single aromatic ethylamine chain (aralylamino group) indicative of a monoamine"

# Example usage (uncomment to test):
# smiles_examples = [
#     "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",  # (S)-dobutamine (True positive)
#     "S(OC1=C(O)C=C([C@@H](O)CNC)C=C1)(O)(=O)=O",  # Epinephrine sulfate (was false negative before)
#     "C1(=C(C=C2C(=C1)C3=C(N2)[C@]4(...",       # A false positive example (likely >500 Da)
# ]
# for smi in smiles_examples:
#     result, reason = is_monoamine(smi)
#     print(f"SMILES: {smi}\nMonoamine: {result}\nReason: {reason}\n")