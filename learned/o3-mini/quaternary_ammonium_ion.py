"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: Quaternary Ammonium Ion

A quaternary ammonium ion is defined as a derivative of ammonium, NH4(+),
in which all four of the hydrogens bonded to nitrogen have been replaced with 
univalent (usually organyl) groups.

Our improved approach is to scan for candidate nitrogen centers that:
  - Have atomic number 7, formal charge +1,
  - Are bonded to exactly 4 atoms and have no attached hydrogens.
  - If the nitrogen is not aromatic, it must be sp3-hybridized.
  - Not all four substituents should be in rings.
  - We then exclude candidates that appear to belong to carnitine-like structures;
    that is, if any carbon neighbor is part of a carboxylate motif (C(=O)[O-]), we skip it.
  - For molecules containing phosphorus (for example, phospholipids) we require that the candidate 
    be within 5 bonds (topological distance) of at least one phosphorus atom.
  - For molecules with no phosphorus and many heavy atoms (>40), we require that the candidate
    have at least one aromatic neighbor.
If any candidate passes all these filters, the function returns True and provides the candidate’s atom index.
Otherwise, it returns False with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a quaternary ammonium ion based on its SMILES string.
    
    The strategy is as follows:
      1. Find candidate nitrogen atoms that:
           - Are nitrogen (atomic num 7) with formal charge +1.
           - Are bonded to exactly 4 other atoms.
           - Have no attached hydrogens (explicit or implicit).
           - If non-aromatic, must be sp3-hybridized.
      2. Discard a candidate if all its substituents are in rings.
      3. Exclude a candidate if any of its carbon substituents is part of a carboxylate group,
         i.e. the carbon is bonded to an oxygen with a double bond and a formal charge -1.
      4. If phosphorus atoms exist in the molecule, ensure the candidate is within 5 bonds of at least one.
      5. In molecules having no phosphorus and being relatively large (>40 atoms),
         require that the candidate has at least one aromatic neighbor.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if a quaternary ammonium ion candidate is found, False otherwise.
       str: A reason message including the candidate’s atom index if found or an explanation if not.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    num_atoms = mol.GetNumAtoms()
    p_atom_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
    
    # Helper function: compute the topological bond distance between two atoms.
    def bond_distance(idx1, idx2):
        try:
            path = Chem.GetShortestPath(mol, idx1, idx2)
            if path:
                return len(path) - 1
        except Exception:
            pass
        return None

    # Helper function: decide if a given candidate should be excluded based on carboxylate-associated carbon.
    def is_attached_carboxylate(candidate):
        # For each neighbor of the candidate nitrogen that is carbon, check if it is part of a carboxylate.
        for nbr in candidate.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                # Look for an oxygen neighbor that is double-bonded with -1 formal charge.
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            if subnbr.GetFormalCharge() == -1:
                                return True
        return False

    # Iterate through all atoms, looking for candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        if atom.GetFormalCharge() != 1:
            continue
        if atom.GetDegree() != 4:
            continue
        if atom.GetTotalNumHs() != 0:
            continue
        # Allow aromatic N candidates even if not sp3; if not aromatic, require sp3-hybridization.
        if not atom.GetIsAromatic():
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
        
        # Check that not all four substituents are in rings.
        neighbors = atom.GetNeighbors()
        if all(nbr.IsInRing() for nbr in neighbors):
            continue

        # Exclude candidate if any one of its carbon substituents appears to be part of a carboxylate.
        if is_attached_carboxylate(atom):
            continue

        candidate_idx = atom.GetIdx()
        passes_extra = False

        # If the molecule contains phosphorus (common in phospholipids),
        # require that the candidate nitrogen be within 5 bonds of at least one phosphorus atom.
        if p_atom_indices:
            for p_idx in p_atom_indices:
                dist = bond_distance(candidate_idx, p_idx)
                if dist is not None and dist <= 5:
                    passes_extra = True
                    break
            if not passes_extra:
                continue
        else:
            # In molecules with no phosphorus and being large (>40 heavy atoms),
            # require that the candidate has at least one aromatic neighbor.
            if num_atoms > 40:
                for nbr in neighbors:
                    if nbr.GetIsAromatic():
                        passes_extra = True
                        break
                if not passes_extra:
                    continue
            else:
                passes_extra = True

        if passes_extra:
            return True, f"Found quaternary ammonium ion at atom index {candidate_idx}"

    return False, "No suitable quaternary ammonium ion found in the molecule"

# Uncomment the following block to run a simple test.
if __name__ == "__main__":
    # Example test: bretylium (should be true)
    test_smiles = "CC[N+](C)(C)Cc1ccccc1Br"
    result, reason = is_quaternary_ammonium_ion(test_smiles)
    print(result, reason)