"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
#!/usr/bin/env python3
"""
Classifies: semisynthetic derivative 
Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis

Improvements over the previous heuristic:
  • Require at least 2 rings, but if the ring system is very complex, lack of chiral centers is less penalized.
  • Adjust the minimal required side-chain modifications (atoms outside the Murcko scaffold) dynamically 
    based on the overall heavy atom count.
  • Enforce a minimum heavy atom count.
Note: This heuristic uses several rules that will not catch all edge cases.
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
    # Parse and check SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure molecule is organic (has carbon atoms) and has enough heavy atoms
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not appear to be organic (no carbon atoms detected)"
    mol_heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    if mol_heavy < 15:
        return False, "Molecule has too few heavy atoms to be a typical natural product derivative"
    
    # Extract Murcko scaffold as a proxy for the natural product core
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_heavy = rdMolDescriptors.CalcNumHeavyAtoms(scaffold)
    side_chain_atoms = mol_heavy - scaffold_heavy
    # Use a dynamic threshold: for smaller molecules require ~10% modification;
    # for larger molecules (>40 heavy atoms) allow a lower ratio (5%)
    modification_threshold = 0.10 if mol_heavy < 40 else 0.05
    side_chain_ratio = side_chain_atoms / mol_heavy if mol_heavy > 0 else 0

    # Count ring systems (overall molecule ring count)
    ring_info = mol.GetRingInfo()
    mol_ring_count = ring_info.NumRings()
    
    # Count chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral = len(chiral_centers)
    
    reasons = []
    # Rule 1. Check for complex ring system.
    if mol_ring_count < 2:
        reasons.append("Insufficient ring complexity (fewer than 2 rings), which is unusual for natural product cores")
    
    # Rule 2. Check stereochemical complexity: if the molecule has only 1 or 2 rings,
    # then require at least one chiral center. For more complex ring systems, allow a low chiral count.
    if mol_ring_count < 3 and num_chiral == 0:
        reasons.append("No chiral centers detected in a low-ring-count system; natural products are usually stereochemically complex")
    
    # Rule 3. Check that a reasonable fraction of atoms lie outside the scaffold.
    if side_chain_ratio < modification_threshold:
        reasons.append("Only minimal modifications detected outside the core scaffold "
                       f"(side chain ratio {side_chain_ratio:.2f} is below the threshold of {modification_threshold:.2f})")
    
    # If any heuristic is violated, label as not semisynthetic derivative.
    if reasons:
        return False, "; ".join(reasons)
    
    # Otherwise, classify as semisynthetic derivative.
    return True, ("Molecule has a sufficiently complex ring system and (if needed) chiral centers, "
                  f"with {side_chain_ratio:.2f} of its heavy atoms in modifications; "
                  "this is consistent with a semisynthetic derivative derived from a natural product core.")

# Example usage:
if __name__ == "__main__":
    # You may test the function with one of the provided SMILES strings.
    test_smiles = "COc1cc2[nH]c3c(C)nccc3c2cc1Br"  # 6-bromoharmine (previously missed due to no chiral centers)
    result, explanation = is_semisynthetic_derivative(test_smiles)
    print(f"Result: {result}\nExplanation: {explanation}")