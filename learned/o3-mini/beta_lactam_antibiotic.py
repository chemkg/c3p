"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam Antibiotic
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
A beta-lactam ring is defined as a four-membered cyclic amide (azetidinone) in which one ring
atom is nitrogen (N) and one of the ring carbons is involved in a C=O (carbonyl) group.
This improved algorithm first “cleans” the input (by removing small counter ions and solvents)
and then searches (via a SMARTS pattern and a manual ring inspection) for a suitable beta-lactam ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    The method:
      1. Parses and sanitizes the molecule.
      2. If multiple fragments are present, picks the largest fragment.
      3. Rejects molecules that contain common metal counter ions.
      4. Searches for a beta-lactam ring using an explicit SMARTS pattern designed to match
         a four-membered cyclic amide (i.e. one nitrogen and a carbonyl-bearing carbon).
      5. If the SMARTS search fails, it loops manually over all 4-membered rings and checks
         that (a) the ring contains exactly one nitrogen atom and (b) one of the carbon atoms in
         the ring carries a C=O bond (with the oxygen not part of the ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a beta-lactam antibiotic, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove salts / waters by taking the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFragments=True)
    if len(frags) > 1:
        # Choose the fragment with the largest number of heavy atoms.
        mol = max(frags, key=lambda m: rdMolDescriptors.CalcHeavyAtomCount(m))
    
    # Reject if explicit metal atoms are found (to avoid counter ions)
    # (common metal symbols found in salts: Na, K, Li, etc.)
    metal_symbols = {"Na", "K", "Li", "Mg", "Ca", "Fe", "Zn"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metal_symbols:
            return False, "Molecule contains metal counter ion(s), likely a salt rather than the active beta-lactam antibiotic."

    # Define a SMARTS pattern for a generic beta-lactam ring.
    # This pattern requires a ring (R) of 4 atoms in which one is nitrogen and
    # one of the carbons has a C=O (carbonyl) attached.
    # Note: The pattern [N;R]1[C;R](=O)[C;R][C;R]1 ~= beta-lactam.
    beta_lactam_smarts = "[N;R]1[C;R](=O)[C;R][C;R]1"
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    if beta_lactam_pattern is None:
        return False, "SMARTS pattern creation failed"

    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, ("Molecule contains a beta-lactam ring (as matched by SMARTS '{}'), "
                      "which is consistent with a beta-lactam antibiotic.".format(beta_lactam_smarts))
    
    # If the SMARTS match did not trigger a hit, do a detailed manual inspection of 4-membered rings.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    for ring in atom_rings:
        if len(ring) == 4:
            n_count = 0
            carbonyl_found = False
            # Loop over the ring atoms
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Count the number of nitrogen atoms in the ring
                if atom.GetAtomicNum() == 7:
                    n_count += 1
                # For carbons, check for a double bond to oxygen outside of the ring (i.e. the carbonyl)
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # Check if bond is a double bond (bond order approximately 2)
                        if bond.GetBondTypeAsDouble() == 2.0:
                            other = bond.GetOtherAtom(atom)
                            # Ensure the other atom is oxygen and not part of the same ring 
                            # (to avoid picking up lactone rings instead of lactams)
                            if other.GetAtomicNum() == 8 and other.GetIdx() not in ring:
                                carbonyl_found = True
                                break
                # (No early exit here; we want to count exactly for this ring.)
            if n_count == 1 and carbonyl_found:
                return True, ("Molecule contains a 4‐membered ring with one nitrogen and a carbonyl carbon, "
                              "which is consistent with a beta–lactam ring.")
    
    return False, "Molecule does not contain a beta-lactam ring."

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with a simple beta-lactam (azetidin-2-one)
    test_smiles = "O=C1CCN1"  
    result, reason = is_beta_lactam_antibiotic(test_smiles)
    print(result, reason)