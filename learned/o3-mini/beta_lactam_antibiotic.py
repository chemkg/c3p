"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam Antibiotic
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
A beta-lactam ring is defined as a four-membered cyclic amide (azetidinone) where one ring atom is nitrogen and 
one of the ring carbons is involved in a carbonyl (C=O) bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    Steps:
      1. Parse the SMILES string. If invalid, return False.
      2. Remove salts and small fragments by selecting the largest fragment.
      3. Reject molecules containing common metal counter ions.
      4. Use a SMARTS pattern to search for a beta-lactam ring.
      5. If SMARTS does not find a match, manually inspect all 4-membered rings for:
         a. Exactly one nitrogen atom.
         b. At least one carbon that has a double-bonded oxygen (a carbonyl) on an exocyclic position.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a beta-lactam antibiotic, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Remove salts and small fragments by taking the largest fragment.
    # Note: The correct parameter is 'sanitizeFrags' not 'sanitizeFragments'.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        # Select fragment with the most heavy atoms.
        mol = max(frags, key=lambda m: rdMolDescriptors.CalcHeavyAtomCount(m))
    
    # Reject molecules that contain common metal counter ions to avoid salts.
    metal_symbols = {"Na", "K", "Li", "Mg", "Ca", "Fe", "Zn"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metal_symbols:
            return False, "Molecule contains metal counter ion(s), likely a salt rather than an active beta-lactam antibiotic."
    
    # Define a SMARTS pattern for the beta-lactam ring:
    # This pattern matches a 4-membered ring: one nitrogen and one C=O edge.
    beta_lactam_smarts = "[N;R]1[C;R](=O)[C;R][C;R]1"
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    if beta_lactam_pattern is None:
        return False, "Failed to create SMARTS pattern for beta-lactam ring."
    
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, ("Molecule contains a beta-lactam ring (matched by SMARTS '{}'), "
                      "consistent with a beta-lactam antibiotic.".format(beta_lactam_smarts))
    
    # If SMARTS didn't match, manually inspect all 4-membered rings in the molecule.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 4:
            n_count = 0
            carbonyl_found = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Count nitrogens in the ring.
                if atom.GetAtomicNum() == 7:
                    n_count += 1
                # For carbons, check if there is a C=O bond with oxygen not part of the ring.
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        if bond.GetBondTypeAsDouble() == 2.0:
                            other = bond.GetOtherAtom(atom)
                            if other.GetAtomicNum() == 8 and other.GetIdx() not in ring:
                                carbonyl_found = True
                                break
            if n_count == 1 and carbonyl_found:
                return True, ("Molecule contains a four-membered ring with one nitrogen and a carbonyl-bearing carbon, "
                              "consistent with a beta-lactam ring.")
    
    return False, "Molecule does not contain a beta-lactam ring."

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with a simple beta-lactam (azetidin-2-one)
    test_smiles = "O=C1CCN1"
    result, reason = is_beta_lactam_antibiotic(test_smiles)
    print(result, reason)