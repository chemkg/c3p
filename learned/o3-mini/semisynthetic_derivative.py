"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
#!/usr/bin/env python3
"""
Classifies: semisynthetic derivative
Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis.
Improved heuristic:
  • Uses Murcko scaffold to define the “natural product core” and requires that it be sufficiently large.
  • Requires at least 3 side‐chain heavy atoms.
  • Uses a dynamic threshold for the fraction of heavy atoms outside the core.
      – For small molecules (<40 heavy atoms): base threshold = 0.10
      – For larger molecules: base threshold = 0.05
  • If no chiral centers are present (which is unusual for natural product cores) the required fraction is increased by 0.05.
  • If halogen(s) are present (commonly introduced by synthesis) the threshold is reduced by 0.02 (with a floor of 0.03).
  • Checks ring complexity with a “bonus” if the side‐chain decoration is very high.
  • Requires that the scaffold (core) should be reasonably complex (>= 6 heavy atoms).
Note: This heuristic still makes mistakes.
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
        str: Explanation/reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule is organic (at least one carbon atom).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not appear to be organic (no carbon atoms detected)"
    
    # Require a minimum number of heavy atoms (non-hydrogen).
    mol_heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    if mol_heavy < 15:
        return False, "Molecule has too few heavy atoms to be a typical natural product derivative"
    
    # Obtain the Murcko scaffold -- a proxy for the natural product core.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_heavy = rdMolDescriptors.CalcNumHeavyAtoms(scaffold)
    # We require that the core be substantial (set arbitrarily at >=6 heavy atoms).
    if scaffold_heavy < 6:
        return False, f"Scaffold is too simple (only {scaffold_heavy} heavy atoms) to represent a natural product core"
    
    # Compute side-chain (decoration) atoms.
    side_chain_atoms = mol_heavy - scaffold_heavy
    side_chain_ratio = side_chain_atoms / mol_heavy if mol_heavy > 0 else 0

    # Dynamic threshold based on overall size.
    base_threshold = 0.10 if mol_heavy < 40 else 0.05

    # Count number of chiral centers.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral = len(chiral_centers)

    # If no chiral centers are present, require a higher fraction modification 
    # (i.e. the natural product core is typically stereochemically complex).
    if num_chiral == 0:
        base_threshold += 0.05

    # Look for halogens – these are common in semisynthetic modifications.
    halogen_numbers = {9, 17, 35, 53}  # F, Cl, Br, I
    has_halogen = any(atom.GetAtomicNum() in halogen_numbers for atom in mol.GetAtoms())
    if has_halogen:
        # Slightly relax the verification of modification when halogen(s) are present.
        base_threshold = max(base_threshold - 0.02, 0.03)
        
    # Count rings.
    mol_ring_count = mol.GetRingInfo().NumRings()
    
    # Initialize list for reasons to reject.
    reasons = []
    
    # Absolute requirement for side-chain atoms: require at least 3 heavy atoms appended.
    if side_chain_atoms < 3:
        reasons.append(f"Only {side_chain_atoms} atoms outside the scaffold; need at least 3 for a semisynthetic derivative")
        
    # Ring complexity rules.
    # Standard expectation: at least 2 rings.
    # However, if there is a single ring then we demand that the side-chain decoration is very high (>=0.30)
    # and that there is some stereochemical complexity.
    if mol_ring_count < 2:
        if not (mol_ring_count == 1 and side_chain_ratio >= 0.30 and num_chiral >= 1):
            reasons.append("Insufficient ring complexity (fewer than 2 rings and no compensating side-chain decoration or stereochemistry)")
    else:
        # For a two‐ring system, having zero chiral centers is unusual for a natural product core.
        if mol_ring_count == 2 and num_chiral == 0:
            reasons.append("A two‐ring system without any chiral centers is unusual for a natural product core")
    
    # Check if the fraction of heavy atoms outside the core meets the threshold.
    if side_chain_ratio < base_threshold:
        reasons.append(f"Modifications outside the core are minimal (side chain ratio {side_chain_ratio:.2f} is below threshold of {base_threshold:.2f})")
    
    # If any criteria are unmet, reject.
    if reasons:
        return False, "; ".join(reasons)
    
    # Otherwise, the molecule is accepted as a (likely) semisynthetic derivative.
    explanation = (
        f"Molecule has {mol_ring_count} ring(s), {num_chiral} chiral center(s), "
        f"{side_chain_atoms} atoms (i.e. {side_chain_ratio:.2f} of heavy atoms) outside the core (threshold: {base_threshold:.2f}), "
        f"and a core of {scaffold_heavy} heavy atoms. This is consistent with a semisynthetic derivative derived from a natural product core."
    )
    return True, explanation

# Example usage: testing with a molecule (you may plug in any of the examples provided)
if __name__ == "__main__":
    # For demonstration, here is the SMILES for TX-1875 (one of the examples)
    test_smiles = r"C(\[C@H]([C@]1([C@H](/C=C/CCC[C@H](C)O)C[C@@H](C1)O)[H])O)=C/C(=O)O"
    result, explanation = is_semisynthetic_derivative(test_smiles)
    print(f"Result: {result}")
    print(f"Explanation: {explanation}")