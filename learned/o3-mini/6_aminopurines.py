"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: 6-aminopurines (any compound having 6-aminopurine (adenine) as part of its structure)
We attempt to detect adenine by:
  1. Looking for the purine ring system (using a SMARTS query).
  2. Checking that one of the ring atoms (a carbon) is substituted (exocyclic) with a primary amine.
  
This two‐step approach is more robust for highly decorated compounds (such as acyl‐CoA derivatives)
that contain adenine as part of their structure.
"""

from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety.
    Our algorithm uses two steps:
      1. Search for a purine ring system using a flexible SMARTS pattern.
      2. Look for a primary exocyclic amine (–NH2) attached to one of the atoms in the purine ring.
    This should capture adenine (or its derivatives) even if it is decorated with other groups.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains the 6-aminopurine (adenine) moiety, False otherwise.
        str: Detailed reason for the classification.
    """
    # Convert SMILES to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Look for the purine ring system.
    # The SMARTS "c1nc2ncnc2n1" is a reasonably flexible pattern for the purine scaffold.
    purine_smarts = "c1nc2ncnc2n1"
    purine_query = Chem.MolFromSmarts(purine_smarts)
    if purine_query is None:
        return False, "Error creating purine substructure query"
    
    purine_matches = mol.GetSubstructMatches(purine_query)
    if not purine_matches:
        return False, "No purine ring system detected"
    
    # Store all atom indices that occur in any purine match.
    purine_atoms = set()
    for match in purine_matches:
        purine_atoms.update(match)
    
    # Step 2. Look for an exocyclic amino group attached to the purine.
    # We iterate over all nitrogen atoms in the molecule which are not in a ring,
    # have exactly one neighbor (i.e. they are terminal) and (by RDKit's counting) have two hydrogens.
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen.
        if atom.GetAtomicNum() != 7:
            continue
        # We require that the amino nitrogen is exocyclic (not part of a ring)
        if atom.IsInRing():
            continue
        # We expect a primary amine: one neighbor and two total hydrogens.
        # (GetTotalNumHs returns the total number of bonded hydrogens, both explicit and implicit.)
        if len(atom.GetNeighbors()) != 1 or atom.GetTotalNumHs() != 2:
            continue
        
        # Get the neighbor atom.
        neighbor = atom.GetNeighbors()[0]
        # Check that the neighbor is aromatic and is part of a purine match.
        if neighbor.GetIsAromatic() and neighbor.GetIdx() in purine_atoms:
            return True, "Molecule contains the 6-aminopurine (adenine) moiety"
    
    return False, "6-aminopurine (adenine) moiety not found in the molecule"

# Example usage (for testing purposes):
# test_smiles = "Cn1cnc(N)c2ncnc12"  # 3-methyladenine (should return True)
# result, reason = is_6_aminopurines(test_smiles)
# print(result, reason)