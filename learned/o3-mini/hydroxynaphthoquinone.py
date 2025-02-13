"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: Hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthoquinone moiety is substituted by at least one hydroxy group.

The approach:
1. Parse the SMILES string.
2. Look for a naphthalene ring system (fused bicyclic aromatic system) using the SMARTS "c1ccc2ccccc2c1".
3. For atoms in the naphthalene core select those that are aromatic carbon atoms carrying a carbonyl – (C(=O)) group.
   We require at least two such carbonyls on the naphthalene core to designate it as a quinone.
4. Check that at least one atom of the naphthalene core carries an –OH substituent.
5. If these conditions are met, return True along with a reason; otherwise, return False with a detailed explanation.

Note: This approach uses substructure searches and is one possible interpretation.
"""

from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    
    A hydroxynaphthoquinone is defined as a compound that has a naphthoquinone core (a naphthalene
    ring system having at least two carbonyl groups as part of the ring) and at least one hydroxy group 
    attached to that core.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a naphthalene core.
    # This pattern looks for two fused benzene rings.
    naphthalene_smarts = "c1ccc2ccccc2c1"
    naphthalene_query = Chem.MolFromSmarts(naphthalene_smarts)
    naph_matches = mol.GetSubstructMatches(naphthalene_query)
    if not naph_matches:
        return False, "No naphthalene (fused aromatic) core found"
    
    # Define SMARTS for an aromatic carbonyl. It matches an aromatic carbon (in a ring)
    # directly double-bonded to oxygen.
    carbonyl_smarts = "[c;R]=[O]"
    carbonyl_query = Chem.MolFromSmarts(carbonyl_smarts)
    
    # Define SMARTS for an aromatic hydroxy substituent, i.e. an –OH attached to an aromatic carbon in a ring.
    hydroxyl_smarts = "[c;R][OH]"
    hydroxyl_query = Chem.MolFromSmarts(hydroxyl_smarts)
    
    # We now iterate over each found naphthalene core (each match is a tuple of atom indices).
    # For each candidate, count how many atoms on the core are part of a carbonyl and how many have a hydroxy group.
    for core in naph_matches:
        core_set = set(core)
        # Count carbonyl groups on the core:
        carbonyl_hits = 0
        carbonyl_matches = mol.GetSubstructMatches(carbonyl_query)
        for match in carbonyl_matches:
            # 'match' is a tuple of indices for the carbonyl pattern.
            # We require that the aromatic C (first atom in the match) belongs to the naphthalene core.
            if match[0] in core_set:
                carbonyl_hits += 1
        if carbonyl_hits < 2:
            # This core does not have the required two carbonyls to be considered a quinone
            continue  # try next naphthalene core if any
        
        # Count hydroxy substituents on the core:
        hydroxyl_hits = 0
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_query)
        for match in hydroxyl_matches:
            # In the pattern, the aromatic atom is first.
            if match[0] in core_set:
                hydroxyl_hits += 1
        
        if hydroxyl_hits < 1:
            # even if we have a quinone core, we need a hydroxy substituent on the naphthalene
            return False, "Found naphthoquinone core (>=2 aromatic carbonyls) but no hydroxy substituent on it"
        
        # If we get here, this naphthalene core qualifies as hydroxynaphthoquinone.
        return True, "Found naphthoquinone core (>=2 carbonyls) with at least one hydroxy substituent on the core"
    
    # If none of the naphthalene cores in the molecule have at least two carbonyl groups, then:
    return False, "No naphthoquinone moiety detected (naphthalene core lacks at least 2 aromatic carbonyl groups)"
    
# (The module can be tested by calling the function with one or more SMILES strings)
# Example: print(is_hydroxynaphthoquinone("OC1=CC(=O)c2ccccc2C1=O"))