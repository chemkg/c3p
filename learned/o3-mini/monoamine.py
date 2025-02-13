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
    A monoamine is defined as an aralylamino compound that contains exactly one amino 
    group connected to an aromatic ring by a two‐carbon chain. In catecholamines the first 
    carbon may carry a hydroxyl substituent. Additional filters:
        - Molecule should not contain a carboxyl group (to avoid amino acid derivatives)
        - Molecular weight should be <= 500 Da.
    
    The algorithm:
      1. Discards molecules with carboxyl groups.
      2. Identifies six-membered aromatic rings and records their carbon atom indices.
      3. For each non‐aromatic nitrogen atom in the molecule, it checks if there is 
         a shortest path of exactly 4 atoms (3 bonds) from an aromatic carbon in a six‐membered ring 
         to that nitrogen, with the two intermediate atoms being non‐aromatic carbons.
      4. Requires exactly one unique nitrogen that forms such an "aralylamino" chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a monoamine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Exclude molecules with a carboxyl group (to avoid amino acid derivatives)
    # This SMARTS pattern finds a carbonyl with an -OH (or deprotonated form)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if carboxyl_pattern and mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Contains carboxyl group, likely an amino acid derivative rather than a monoamine"
    
    # Filter by molecular weight: monoamines usually are small (<= 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} Da is too high for a typical monoamine"
    
    # Pre-calculate ring information: record indices of carbons that are in a six-membered ring.
    ring_info = mol.GetRingInfo()
    six_membered_aromatic_atoms = set()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            # Check if the ring is aromatic (most monoamine precursors have benzene rings)
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                # Record only carbon atoms
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() == "C":
                        six_membered_aromatic_atoms.add(idx)
    
    if not six_membered_aromatic_atoms:
        return False, "No six-membered aromatic ring found"
    
    # Now search for an 'aralylamino' chain.
    # We require a path with four atoms: aromatic carbon, linker1, linker2, and an amino nitrogen.
    valid_nitrogen_matches = set()
    num_atoms = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        # Consider only nitrogen atoms that are not aromatic (to avoid ring N's)
        if atom.GetSymbol() != "N" or atom.GetIsAromatic():
            continue
        n_idx = atom.GetIdx()
        # For each candidate six-membered aromatic carbon
        for aro_idx in six_membered_aromatic_atoms:
            # Get the shortest path (as a tuple of atom indices) from the aromatic carbon to this nitrogen.
            # If no path exists, skip.
            path = Chem.rdmolops.GetShortestPath(mol, aro_idx, n_idx)
            # We require exactly 3 bonds i.e. 4 atoms in the path.
            if len(path) != 4:
                continue
            # Validate the chain pattern:
            # - The starting atom must be aromatic (and is in our six_membered_aromatic_atoms by definition)
            # - The two linker atoms (path[1] and path[2]) must be carbons and should not be aromatic.
            linker1 = mol.GetAtomWithIdx(path[1])
            linker2 = mol.GetAtomWithIdx(path[2])
            if linker1.GetSymbol() != "C" or linker2.GetSymbol() != "C":
                continue
            if linker1.GetIsAromatic() or linker2.GetIsAromatic():
                continue
            # Accept the chain if all checks pass.
            valid_nitrogen_matches.add(n_idx)
    
    # Require exactly one unique nitrogen involved in a valid chain.
    if len(valid_nitrogen_matches) == 0:
        return False, "No valid aromatic ethylamine chain found"
    elif len(valid_nitrogen_matches) > 1:
        return False, f"Found {len(valid_nitrogen_matches)} distinct aromatic ethylamine chain(s); expected exactly one"
    
    return True, "Found a single aromatic ethylamine chain (aralylamino group) indicative of a monoamine"

# Example usage (uncomment the lines below to test):
# test_smiles = [
#     "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",  # (S)-dobutamine (should be True)
#     "[C@@H]([C@@H](N)C)(O)C1=CC(O)=C(C=C1)O",      # alpha-methylnoradrenaline (should be True now)
#     "C(CNCCCCCCNCCC1=CC=CC=C1)C2=CC(O)=C(C=C2)O",   # dopexamine (has 2 chains -> False)
#     "N[C@@H](CSc1cc(C[C@H](N)C(O)=O)cc(O)c1O)C(O)=O",  # Cysteinyldopa (should be False)
# ]
# for smi in test_smiles:
#     result, reason = is_monoamine(smi)
#     print(f"SMILES: {smi}\nMonoamine: {result}\nReason: {reason}\n")