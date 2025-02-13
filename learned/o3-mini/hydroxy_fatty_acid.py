"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: Hydroxy fatty acid
Definition: A fatty acid is defined as a molecule possessing a carboxylic acid (or its deprotonated form) 
attached to a predominantly aliphatic chain and carrying one or more hydroxyl substituents (not part of the carboxylic acid).
This implementation uses additional heuristics (e.g. absence of nitrogen, no rings, sufficiently high carbon content)
to reduce false positives.
"""

from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.

    A hydroxy fatty acid is defined as a fatty acid with a carboxyl group (or carboxylate)
    on a predominantly aliphatic (acyclic, mainly carbon) chain that additionally carries one 
    or more free hydroxyl (OH) substituents that are not part of the carboxyl moiety.

    Additional heuristics applied:
      - The molecule must be acyclic (to avoid glycosides, peptides, and other ring‐containing molecules).
      - The molecule should have no nitrogen atoms.
      - The ratio of carbon atoms to heavy atoms must be high (>=0.7) to ensure a predominantly aliphatic character.
      - There must be at least 5 carbon atoms.
      - The carboxyl group is recognized both in its acid “[CX3](=O)[OX2H1]” and deprotonated “[CX3](=O)[O-]” forms.
      - At least one hydroxyl group (SMARTS "[OX2H]") must be present that is not part of a carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a hydroxy fatty acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Reject molecules containing nitrogen since fatty acids typically do not contain N.
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Nitrogen present; likely not a fatty acid."

    # Reject molecules containing rings (most fatty acids are acyclic linear chains).
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; likely not a linear fatty acid."

    # Define SMARTS patterns for a carboxyl group in both acid and deprotonated forms.
    carboxyl_acid_smarts = "[CX3](=O)[OX2H1]"
    carboxylate_smarts   = "[CX3](=O)[O-]"
    carboxyl_acid = Chem.MolFromSmarts(carboxyl_acid_smarts)
    carboxylate   = Chem.MolFromSmarts(carboxylate_smarts)
    
    # Check if the molecule contains either a carboxylic acid or a carboxylate group.
    if not (mol.HasSubstructMatch(carboxyl_acid) or mol.HasSubstructMatch(carboxylate)):
        return False, "No carboxyl group found; not a fatty acid."

    # Gather the oxygen atom indices that are part of the carboxyl function.
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_acid) + mol.GetSubstructMatches(carboxylate)
    carboxyl_oxygens = set()
    # In our SMARTS the oxygen is the second atom (index 1).
    for match in carboxyl_matches:
        if len(match) > 1:
            carboxyl_oxygens.add(match[1])

    # Define the SMARTS for a free hydroxyl group.
    hydroxyl_smarts = "[OX2H]"
    hydroxyl_query = Chem.MolFromSmarts(hydroxyl_smarts)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_query)
    
    # Check if at least one hydroxyl group is present that is not part of a carboxyl group.
    free_hydroxy_found = False
    for match in hydroxyl_matches:
        oh_idx = match[0]
        if oh_idx not in carboxyl_oxygens:
            free_hydroxy_found = True
            break
    if not free_hydroxy_found:
        return False, "Fatty acid found but no additional free hydroxyl substituent detected."

    # Additional heuristic: ensure the molecule is predominantly aliphatic.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if carbon_count < 5:
        return False, "Not enough carbon atoms to be a fatty acid."
    if heavy_atom_count > 0 and (carbon_count / heavy_atom_count) < 0.7:
        return False, "Molecule is not predominantly aliphatic."

    return True, ("Contains a carboxyl (or carboxylate) group on an aliphatic chain and a free hydroxyl substituent, "
                  "classifying it as a hydroxy fatty acid.")

# Example usage (uncomment the following line to test):
# print(is_hydroxy_fatty_acid("C[C@H](O)CCCCC(O)=O"))  # Expected: (6S)-6-hydroxyheptanoic acid