"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: Anilide
Definition: Anilides are aromatic amides formed by acylation of aniline.
That is, the amide group (C(=O)N) must have the acyl carbon
attached to the nitrogen while the nitrogen is directly bonded to exactly one aromatic (benzene) carbon.
Optionally, it may have one additional heavy (non-hydrogen) aliphatic substituent.
This code improves previous performance by stricter counting of neighbors.
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
            # Verify that all atoms in the ring are carbon and aromatic.
            if all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 and mol.GetAtomWithIdx(i).GetIsAromatic() 
                   for i in ring):
                return True
    return False

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    Anilides are aromatic amides formed by the acylation of aniline.
    In our criteria the amide nitrogen must be bonded to:
      - the acyl carbon (of the C(=O) group) and
      - exactly one aromatic (benzene) substituent,
      - optionally one additional heavy (non-hydrogen) aliphatic substituent.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an anilide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Define amide pattern:
    # Here the SMARTS ensures that we have a carbon with a double-bonded oxygen and a single bond to a nitrogen.
    amide_pattern = Chem.MolFromSmarts("[#6](=O)-[#7]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found"
    
    # Evaluate each amide candidate
    for match in amide_matches:
        # our SMARTS returns (acyl carbon, amide nitrogen)
        acyl_c_idx, nitrogen_idx = match[0], match[1]
        n_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Get all heavy (non-hydrogen) neighbors of the nitrogen.
        heavy_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # Exclude the acyl carbon (which is part of the amide bond)
        other_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetIdx() != acyl_c_idx]
        
        # Count how many of these neighbors are aromatic and in a benzene ring.
        aromatic_count = 0
        for nbr in other_neighbors:
            if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic() and is_in_benzene(nbr, ring_info, mol):
                aromatic_count += 1
        
        # The allowed heavy neighbor count (aside from acyl carbon) should be either 1 or 2.
        # If there is only 1 heavy neighbor then it must be aromatic (typical for acetanilide).
        # If there are 2 heavy neighbors then one should be aromatic and the other aliphatic.
        num_other = len(other_neighbors)
        if (num_other == 1 and aromatic_count == 1) or (num_other == 2 and aromatic_count == 1):
            return True, ("Found an amide bond (C(=O)N) where the nitrogen is bound to exactly one aromatic benzene "
                          "ring substituent and acceptable additional substituents, indicative of an anilide structure.")
    
    # If none of the amide candidates satisfy the criteria then it is not an anilide.
    return False, "No amide bond with the nitrogen bound to a single benzene substituent (with allowed extra groups) was found."

# Example testing (commented out):
# test_smiles_list = [
#     "CCCOCCN(C(=O)CCl)c1c(CC)cccc1CC",  # Pretilachlor (true positive expected)
#     "CCC(=O)N(C1CCN(CC1)CCc1ccccc1)c1ccccc1",  # fentanyl (should be false, as N has two aromatic groups)
# ]
#
# for s in test_smiles_list:
#     result, reason = is_anilide(s)
#     print(s, "->", result, ";", reason)