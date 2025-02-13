"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
#!/usr/bin/env python3
"""
Classifies: semisynthetic derivative (any organic molecular entity derived from a natural product by partial chemical synthesis)
Heuristic:
  1. The molecule must be valid and contain carbon.
  2. A Murcko scaffold is extracted; natural products (and their semisynthetic derivatives) tend to have complex ring systems.
  3. Most natural products are stereochemically rich so the presence of chiral centers is expected.
  4. Semisynthetic derivatives are typically modified relative to the natural core. A significant number of atoms
     outside of the scaffold (e.g. >10% of the total heavy atoms) is taken as tentative evidence of a modification.
Note that this function uses several heuristic rules that will not catch all cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is (likely) a semisynthetic derivative, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule is organic (has carbon atoms)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not appear to be organic (no carbon atoms detected)"
    
    # Extract the Murcko scaffold; natural products often have complex polycyclic cores.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    
    # Calculate the number of heavy atoms (non-hydrogen) in full molecule and its scaffold.
    mol_heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    scaffold_heavy = rdMolDescriptors.CalcNumHeavyAtoms(scaffold)
    side_chain_atoms = mol_heavy - scaffold_heavy
    side_chain_ratio = side_chain_atoms / mol_heavy if mol_heavy > 0 else 0
    
    # Count rings in the full molecule (complex cores typically have several rings)
    mol_ring_count = mol.GetRingInfo().NumRings()
    
    # Count chiral centers (found via '@' symbols in SMILES); most natural products are stereochemically complex.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral = len(chiral_centers)
    
    # Heuristic checks:
    # 1. The molecule should have a significant ring system (at least 2 rings is a minimal requirement).
    # 2. It should have one or more chiral centers.
    # 3. There should be a significant amount of substitution/modification relative to the core;
    #    here we check if at least 10% of heavy atoms lie outside the scaffold.
    reasons = []
    if mol_ring_count < 2:
        reasons.append("Insufficient ring complexity (fewer than 2 rings), which is unusual for natural product cores")
    if num_chiral == 0:
        reasons.append("No chiral centers detected, while natural products are typically stereochemically complex")
    if side_chain_ratio < 0.10:
        reasons.append("Only minimal modifications detected outside the core scaffold (side chain ratio < 10%)")
    
    # If any heuristic is not met, we assume it is not a semisynthetic derivative.
    if reasons:
        return False, "; ".join(reasons)
    
    # If all checks are met, we tentatively classify the molecule as a semisynthetic derivative.
    return True, ("Molecule has a complex polycyclic scaffold with stereocenters and a significant proportion of extra atoms "
                  "indicating modifications to a natural product core; this is consistent with a semisynthetic derivative.")