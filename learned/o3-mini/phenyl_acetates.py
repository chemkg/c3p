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
     six-membered ring that is an isolated benzene ring: all atoms in that ring are aromatic 
     carbons and none of the ring atoms is shared with any other ring.
  3. Reject molecules that are likely too complex (here, with molecular weight >350 Da) to be 
     considered simple phenyl acetates.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester (CH3C(=O)-O-) where the ester oxygen is attached
    directly to an aromatic carbon, and that aromatic carbon must belong to an isolated six-membered
    benzene ring (i.e. no fusion with other rings). Additionally, the overall molecular weight must be
    in the range for a relatively simple phenyl acetate (<=350 Da).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a simple phenyl acetate, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute molecular weight; if too heavy, consider the molecule too complex.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.2f} too high for a simple phenyl acetate"
    
    # Define a SMARTS pattern for the acetate ester group:
    # The pattern enforces:
    # - an aromatic carbon (c) bound to an oxygen (O)
    # - that oxygen bound to a carbonyl (C=O) which is bound to a CH3.
    acetate_pattern = Chem.MolFromSmarts("cOC(=O)[CH3]")
    if acetate_pattern is None:
        return False, "Error processing SMARTS pattern"
    
    # Find all substructure matches for the acetate group in the molecule.
    matches = mol.GetSubstructMatches(acetate_pattern)
    if not matches:
        return False, "No phenyl acetate substructure found"
    
    # Get all ring information (each ring is a tuple of atom indices)
    rings = mol.GetRingInfo().AtomRings()
    
    # For each match, check the aromatic carbon (the first atom in the match) is part of an isolated benzene ring.
    for match in matches:
        aromatic_atom_idx = match[0]  # the aromatic carbon attached to oxygen
        # Iterate over rings that contain this aromatic carbon.
        for ring in rings:
            if aromatic_atom_idx in ring and len(ring) == 6:
                # First, verify that every atom in this ring is aromatic and is a carbon.
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    # Next, check that none of these ring atoms is part of any other ring (i.e. isolated ring).
                    isolated_ring = True
                    for idx in ring:
                        count = sum(1 for r in rings if idx in r)
                        if count > 1:
                            isolated_ring = False
                            break
                    if isolated_ring:
                        return True, "Found phenyl acetate group (acetate ester directly attached to an isolated benzene ring)"
    
    return False, "Phenyl acetate substructure found, but not attached to an isolated benzene ring"