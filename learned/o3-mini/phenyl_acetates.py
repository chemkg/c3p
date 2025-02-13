"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid 
with the hydroxy group of any phenol.
Improvement strategy:
  1. Use a SMARTS pattern that exactly requires the acetate ester part as CH3C(=O)–.
  2. Check that the oxygen is directly attached to an aromatic carbon that belongs to at least one
     six-membered ring which is a simple benzene ring (i.e. all atoms in that ring are aromatic carbons).
  3. Reject molecules that are likely too complex (here, with molecular weight >350 Da) to be considered
     simple phenyl acetates.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester (CH3C(=O)–O–) where the ester oxygen is attached
    directly to a phenol-like (aromatic) ring. The method uses a strict SMARTS to find a CH3C(=O) group
    connected by an oxygen to an aromatic carbon, then checks that this aromatic carbon belongs to 
    at least one 6-membered benzene ring. Additionally, the overall weight must be within range for a 
    simple phenyl acetate (here set as <=350 Da).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phenyl acetate, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the exact molecular weight; reject if too heavy for a "simple" phenyl acetate.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.2f} too high for a simple phenyl acetate"
    
    # Define a SMARTS pattern for a phenyl acetate group: the pattern enforces:
    # - an aromatic carbon (c) bound to an oxygen (O)
    # - that oxygen is bound to a carbonyl carbon (C(=O)) which in turn is bound to a CH3 group.
    acetate_pattern = Chem.MolFromSmarts("cOC(=O)[CH3]")
    if acetate_pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Find all substructure matches for the acetate pattern.
    matches = mol.GetSubstructMatches(acetate_pattern)
    if not matches:
        return False, "No phenyl acetate substructure found"

    # Get ring information for checking the aromatic benzene ring
    ring_info = mol.GetRingInfo().AtomRings()

    # For each match, verify that the aromatic carbon (first atom in the match) is part of a 6-membered ring
    # composed entirely of aromatic carbons.
    for match in matches:
        aromatic_atom_idx = match[0]  # this should be the aromatic carbon attached to the oxygen
        # Check over all rings that contain this atom
        for ring in ring_info:
            if aromatic_atom_idx in ring and len(ring) == 6:
                # Verify that all atoms in the ring are aromatic carbons (atomic number 6)
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() 
                       for idx in ring):
                    return True, "Found phenyl acetate group (acetate ester directly attached to a benzene ring)"
    
    # If none of the acetate pattern matches are connected to a simple benzene ring, reject it.
    return False, "Phenyl acetate substructure found, but not attached to a simple benzene ring"