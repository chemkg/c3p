"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: Anilide
Definition: Anilide is any aromatic amide obtained by acylation of aniline.
That is, the amide (C(=O)N) bond must have the acyl carbon attached to the nitrogen and
the nitrogen must be directly bonded to an aromatic (benzene-like) carbon.
Additionally, aside from the carbonyl carbon, the nitrogen should have exactly one aromatic substituent,
allowing at most one extra (aliphatic) heavy substituent.
"""
from rdkit import Chem

def is_in_benzene(atom, ring_info, mol):
    """
    Checks if the given atom is part of a benzene ring.
    A benzene ring is defined as a six-membered ring in which every atom is carbon and aromatic.
    
    Args:
        atom: an RDKit Atom object.
        ring_info: list of rings (each a tuple of atom indices) from mol.GetRingInfo().
        mol: the RDKit Mol object.
        
    Returns:
        True if the atom is in a benzene ring, False otherwise.
    """
    idx = atom.GetIdx()
    for ring in ring_info:
        if idx in ring and len(ring) == 6:
            # Check that all atoms in the ring are carbons and aromatic.
            if all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 and mol.GetAtomWithIdx(i).GetIsAromatic() 
                   for i in ring):
                return True
    return False

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    Anilides are aromatic amides formed by acylation of aniline,
    meaning that the amide nitrogen (N in C(=O)N) must be directly bonded to a benzene ring.
    To be more selective, we require that aside from the acyl carbon, the nitrogen has exactly one
    aromatic (benzene) substituent and at most one additional heavy neighbor (assumed aliphatic).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Define an amide bond pattern: a carbon bonded to a carbonyl oxygen (=O) and to a nitrogen.
    # This simple SMARTS "[#6](=O)[#7]" will catch most C(=O)N fragments.
    amide_pattern = Chem.MolFromSmarts("[#6](=O)[#7]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found"
    
    # Evaluate each found amide bond candidate
    for match in amide_matches:
        # By our SMARTS, the match returns: (acyl carbon, oxygen, amide nitrogen)
        acyl_c_idx = match[0]
        nitrogen_idx = match[2]
        n_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Get all heavy (non-hydrogen) neighbors of the nitrogen.
        heavy_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # Exclude the acyl carbon (which is part of the amide bond)
        other_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetIdx() != acyl_c_idx]
        
        # Count how many of the remaining heavy neighbors are carbons in a benzene ring.
        benzene_count = 0
        for nbr in other_neighbors:
            if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic() and is_in_benzene(nbr, ring_info, mol):
                benzene_count += 1
        
        # For true anilides (from acylation of aniline), we expect the following:
        #   - There must be exactly one aromatic (benzene) substituent.
        #   - In total, aside from the acyl carbon, at most one extra heavy substituent is acceptable.
        if benzene_count == 1 and len(other_neighbors) <= 2:
            return (True, "Found an amide bond with the N attached to a benzene ring (6-membered aromatic) "
                          "indicative of anilide structure.")
    
    return (False, "No amide bond with the nitrogen bound to a benzene ring (with appropriate substitution) was found.")

# For testing purposes you can run:
# test_smiles_list = [
#    "C1C2CC3(CC1CC(C3)(C2)C(=O)OCC(=O)NC4=CC=CC=C4OC(F)F)",  # a true anilide example with additional groups
#    "CCCOCCN(C(=O)CCl)c1c(CC)cccc1CC",  # Pretilachlor (true positive expected)
#    "CCC(=O)N(C1CCN(CC1)CCc1ccccc1)c1ccccc1",  # fentanyl (should be false, as N has two aromatic groups)
# ]
# for s in test_smiles_list:
#    print(s, "->", is_anilide(s))