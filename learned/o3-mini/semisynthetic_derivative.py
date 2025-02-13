"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
#!/usr/bin/env python3
"""
Classifies: semisynthetic derivative
Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis.
Revision: adjusted thresholds to address a low F1 score.
Improvements include:
  • A dynamic threshold for the fraction of heavy atoms outside the core Murcko scaffold (side chains).
  • An absolute requirement on the number of side‐chain atoms.
  • A relaxed ring count if the side‐chain decorations are very high.
  • A “bonus” if halogen(s) are present (commonly introduced by chemical synthesis).
Note: No heuristic will work perfectly; this is only an approximate approach.
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
        bool: True if the molecule is (likely) a semisynthetic derivative, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Must be an organic molecule (i.e. contain at least one carbon atom).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not appear to be organic (no carbon atoms detected)"
    
    # Ensure the molecule has a substantial number of heavy atoms.
    mol_heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    if mol_heavy < 15:
        return False, "Molecule has too few heavy atoms to be a typical natural product derivative"
    
    # Get the Murcko scaffold as a proxy for the natural product core.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_heavy = rdMolDescriptors.CalcNumHeavyAtoms(scaffold)
    side_chain_atoms = mol_heavy - scaffold_heavy
    # Calculate the ratio of atoms outside the scaffold.
    side_chain_ratio = side_chain_atoms / mol_heavy if mol_heavy > 0 else 0

    # Define a base dynamic threshold:
    # For molecules with less than 40 heavy atoms require ~10% modification,
    # for larger molecules allow roughly 5%.
    modification_threshold = 0.10 if mol_heavy < 40 else 0.05

    # Look for halogen atoms; their presence is common in semisynthesis.
    halogen_numbers = {9, 17, 35, 53}  # F, Cl, Br, I
    has_halogen = any(atom.GetAtomicNum() in halogen_numbers for atom in mol.GetAtoms())
    # If halogen is present, modest modifications might be enough.
    if has_halogen:
        modification_threshold = max(modification_threshold - 0.02, 0.03)

    # Count overall rings in the molecule.
    ring_info = mol.GetRingInfo()
    mol_ring_count = ring_info.NumRings()

    # Count chiral centers (including those not assigned, to capture stereochemical complexity).
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral = len(chiral_centers)

    # Initialize a list to accumulate reasons for rejection.
    reasons = []
    
    # Absolute requirement: require that at least 4 atoms (side chains) are appended.
    if side_chain_atoms < 4:
        reasons.append(f"Only {side_chain_atoms} atoms outside the scaffold; need at least 4 for a semisynthetic derivative")
    
    # Rule for ring complexity:
    # Normally expect at least 2 rings. However, if the side-chain decoration is high (>=0.30)
    # and there is at least one chiral center, then 1-ring molecules might be semisynthetic derivatives.
    if mol_ring_count < 2:
        if side_chain_ratio < 0.30 or num_chiral < 1:
            reasons.append("Insufficient ring complexity (fewer than 2 rings) and inadequate side-chain decoration or stereochemistry")
    else:
        # For molecules with only 2 rings, require at least one chiral center for added stereochemical complexity.
        if mol_ring_count == 2 and num_chiral == 0:
            reasons.append("A two-ring system without any chiral centers is unusual for a natural product core")
    
    # Check that enough of the heavy atoms lie outside the core.
    if side_chain_ratio < modification_threshold:
        reasons.append(f"Modifications outside the core are minimal (side chain ratio {side_chain_ratio:.2f} is below threshold of {modification_threshold:.2f})")
    
    # If any of the heuristics are unmet, do not classify as semisynthetic derivative.
    if reasons:
        return False, "; ".join(reasons)
    
    # Otherwise, the molecule passes our semisynthetic derivative criteria.
    return True, (
        f"Molecule has {mol_ring_count} ring(s), {num_chiral} chiral center(s), "
        f"and {side_chain_ratio:.2f} of its heavy atoms outside the core (with a threshold of {modification_threshold:.2f}); "
        "this is consistent with a semisynthetic derivative derived from a natural product core."
    )

# Example usage
if __name__ == "__main__":
    # You can test with one of the examples such as TX-1875:
    test_smiles = r"C(\[C@H]([C@]1([C@H](/C=C/CCC[C@H](C)O)C[C@@H](C1)O)[H])O)=C/C(=O)O"  # TX-1875
    result, explanation = is_semisynthetic_derivative(test_smiles)
    print(f"Result: {result}\nExplanation: {explanation}")