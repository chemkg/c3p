"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue – a molecule that can substitute for a normal nucleobase in nucleic acids.
Heuristic criteria used:
  1. Valid SMILES and molecular weight in the 80–350 Da range.
  2. Contains at least two nitrogen atoms.
  3. Contains at least one aromatic ring that includes a nitrogen atom.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
    
    Heuristic criteria:
      - Must be a valid molecule.
      - Molecular weight between 80 and 350 Da.
      - Contains at least 2 nitrogen atoms.
      - Contains at least one aromatic ring which has at least one nitrogen atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fits the criteria, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.1f} Da is not in the typical range for nucleobase analogues (80-350 Da)"
    
    # Count total nitrogen atoms
    num_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if num_nitrogen < 2:
        return False, f"Found only {num_nitrogen} nitrogen atom(s); nucleobase analogues typically have 2 or more"
    
    # Retrieve ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "Molecule has no rings; nucleobase analogues are heterocyclic"
    
    # Check for at least one aromatic ring that contains a nitrogen atom
    found_aromatic_heterocycle = False
    for ring in atom_rings:
        # Check if all atoms in the ring are aromatic
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Count nitrogen atoms in this ring
            n_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_in_ring >= 1:
                found_aromatic_heterocycle = True
                break

    if not found_aromatic_heterocycle:
        return False, "Molecule does not have an aromatic heterocycle containing nitrogen"
    
    return True, "Molecule fits heuristic criteria for a nucleobase analogue"