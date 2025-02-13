"""
Classifies: CHEBI:33856 aromatic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: Aromatic Amino Acid
Definition: An amino acid whose structure includes an aromatic ring.
Improved criteria:
  - Must contain a carboxyl (or carboxylate) group in the amino acid backbone.
  - Must contain an amino group (which may be N‐substituted) attached at the α–carbon.
  - Must contain at least one aromatic ring.
  - Must be a small molecule (heavy atom count <= 30) to avoid peptides and larger molecules.
The approach uses a backbone SMARTS that looks for the fragment:
   [NX3;!$(NC(=O))][C;H1](C(=O)[O;H1,-])
which captures the connectivity of most (non–glycine) amino acids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid.
    
    To be classified as an aromatic amino acid, the molecule must:
      - Have an amino acid backbone with a nitrogen directly bound to an α–carbon
        that carries a carboxyl (or carboxylate) group.
      - Contain at least one aromatic ring.
      - Not be too large (heavy atom count <= 30).
    
    This method uses a backbone SMARTS pattern:
      [NX3;!$(NC(=O))][C;H1](C(=O)[O;H1,-])
    which will match both free and N–substituted amino groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets criteria for an aromatic amino acid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are correct.
    mol = Chem.AddHs(mol)
    
    # Check molecular size to avoid peptides and larger biomolecules.
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms > 30:
        return False, f"Molecule is too large (has {heavy_atoms} heavy atoms) to be a single amino acid"
    
    # Check for at least one aromatic ring.
    ring_info = mol.GetRingInfo()
    aromatic_ring_found = False
    for ring in ring_info.AtomRings():
        if len(ring) < 5:
            continue
        # Verify every atom in the ring is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_found = True
            break
    if not aromatic_ring_found:
        return False, "No aromatic ring found"
    
    # SMARTS pattern to capture the amino acid backbone:
    # [NX3;!$(NC(=O))] : A trivalent nitrogen not attached to a carbonyl (thus not an amide).
    # [C;H1]         : An α–carbon carrying exactly one hydrogen (typical for non–glycine residues).
    # (C(=O)[O;H1,-]) : A carboxyl group (which can be protonated or deprotonated).
    amino_acid_smarts = "[NX3;!$(NC(=O))][C;H1](C(=O)[O;H1,-])"
    aa_pattern = Chem.MolFromSmarts(amino_acid_smarts)
    if aa_pattern is None:
        return False, "Error in parsing amino acid SMARTS pattern"
    
    if not mol.HasSubstructMatch(aa_pattern):
        return False, "No amino acid backbone (N–α–C–(C=O)[O]) fragment found"
    
    return True, "Molecule contains an amino acid backbone and an aromatic ring, with appropriate size."

# Example usage (this example can be removed or commented out if using as a module)
if __name__ == "__main__":
    # Test examples: You can try with N-methyl-D-dopa and D-phenylalanine
    examples = {
        "N-methyl-D-dopa": "CN[C@H](CC1=CC=C(O)C(O)=C1)C(O)=O",
        "D-phenylalanine": "N[C@H](Cc1ccccc1)C(O)=O"
    }
    for name, smi in examples.items():
        result, reason = is_aromatic_amino_acid(smi)
        print(f"{name}: {result}\n  Reason: {reason}\n")