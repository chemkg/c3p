"""
Classifies: CHEBI:33856 aromatic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: Aromatic Amino Acid 
Definition: An amino acid whose structure includes an aromatic ring.
Improved criteria:
  - Must contain a carboxylic acid (or carboxylate) group.
  - Must contain at least one amino (N–H) group.
  - Must contain at least one aromatic ring.
  - Molecular size (heavy atom count) is checked to avoid peptides.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid.
    
    To be classified as an aromatic amino acid, the molecule must:
      - Have a carboxyl group (protonated or deprotonated).
      - Have at least one amino group (as evidenced by a nitrogen with at least one attached hydrogen).
      - Contain at least one aromatic ring.
      - Not be too large (which would indicate a peptide or larger compound).
    
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
    
    # Add explicit hydrogens so that hydrogen counts (especially on nitrogen) are correct.
    mol = Chem.AddHs(mol)
    
    # Check molecular size to avoid peptides and larger biomolecules.
    # Typical free amino acids have relatively small heavy atom counts; here we use 30 as an upper bound.
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms > 30:
        return False, f"Molecule is too large (has {heavy_atoms} heavy atoms) to be a single amino acid"
    
    # Check for a carboxyl group.
    # Two SMARTS patterns: one for the protonated acid and one for the carboxylate.
    carboxyl_acid_smarts = "[CX3](=O)[OX2H]"
    carboxylate_smarts   = "[CX3](=O)[O-]"
    carboxyl_acid = Chem.MolFromSmarts(carboxyl_acid_smarts)
    carboxylate   = Chem.MolFromSmarts(carboxylate_smarts)
    
    has_carboxyl = mol.HasSubstructMatch(carboxyl_acid) or mol.HasSubstructMatch(carboxylate)
    if not has_carboxyl:
        return False, "No carboxylic acid (or carboxylate) group found"
    
    # Check for an amino group.
    # Instead of relying on a generic SMARTS that may miss implicit hydrogens,
    # we iterate over nitrogen atoms and check that at least one has a hydrogen attached.
    amino_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # nitrogen
            # Get total number of attached hydrogens (explicit + implicit since we added Hs already)
            if atom.GetTotalNumHs() >= 1:
                # To avoid detecting an amide nitrogen (which usually is directly bound to a carbonyl),
                # we check if any neighbor is a carbon with a double bond to oxygen.
                is_amide = False
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6:
                        for bond in mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()):
                            pass  # placeholder if bond info was needed
                        # Look for a C=O (carbonyl) on that carbon.
                        for cnbr in nbr.GetNeighbors():
                            if cnbr.GetAtomicNum() == 8:
                                bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), cnbr.GetIdx())
                                if bond and bond.GetBondTypeAsDouble() == 2.0:
                                    is_amide = True
                                    break
                        if is_amide:
                            break
                # We only accept the nitrogen if it is not in an amide environment.
                if not is_amide:
                    amino_found = True
                    break
    if not amino_found:
        return False, "No free amino group (or N-H moiety) found"
    
    # Check for at least one aromatic ring.
    ring_info = mol.GetRingInfo()
    aromatic_ring_found = False
    for ring in ring_info.AtomRings():
        # Require a minimum size; many meaningful rings are at least 5 atoms.
        if len(ring) < 5:
            continue
        # Verify every atom in the ring is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_found = True
            break
    if not aromatic_ring_found:
        return False, "No aromatic ring found"
    
    return True, "Molecule contains a carboxyl group, a free amino group, and an aromatic ring."

# Example usage (can be removed or commented out if used as a module):
if __name__ == "__main__":
    # Test with D-phenylalanine (should be classified as aromatic amino acid)
    test_smiles = "N[C@H](Cc1ccccc1)C(O)=O"
    result, reason = is_aromatic_amino_acid(test_smiles)
    print(result, reason)