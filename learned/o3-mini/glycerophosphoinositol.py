"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: glycerophosphoinositol
Definition: Any glycerophospholipid having the polar alcohol inositol 
esterified to the phosphate group at the sn-3 position of the glycerol backbone.
Strategy:
  – Instead of separately checking for a glycerol backbone, a phosphate and an inositol ring,
    we require that at least one phosphate atom is directly bridging an oxygen that belongs to
    a glycerol substructure and another oxygen that belongs to an inositol substructure.
  – This is done by (a) gathering all atoms that are parts of a glycerol backbone (via a
    simplified SMARTS "C(O)C(O)CO") and similarly for a myo‐inositol ring (via SMARTS "C1C(O)C(O)C(O)C(O)C1O")
    and then (b) iterating over all phosphate fragments ([P](=O)(O)(O)) in the molecule. For each phosphate
    we check its attached oxygens: one must “lead into” a glycerol fragment and another into an inositol fragment.
"""

from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    Here we require that the molecule contains a phosphate group that bridges a glycerol moiety (sn-3)
    with an inositol ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a glycerophosphoinositol.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the inositol ring and glycerol backbone.
    # We use simplified patterns.
    inositol_smarts = "C1C(O)C(O)C(O)C(O)C1O"   # myo-inositol (6-membered ring with OH groups)
    glycerol_smarts = "C(O)C(O)CO"                # simplified glycerol backbone (two secondary alcohols and one primary)
    
    inositol_pat = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pat is None:
        return False, "Error in inositol SMARTS pattern"
    
    glycerol_pat = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol_pat is None:
        return False, "Error in glycerol SMARTS pattern"
    
    # Find all substructure matches for the glycerol and inositol fragments.
    glycerol_matches = mol.GetSubstructMatches(glycerol_pat)
    inositol_matches = mol.GetSubstructMatches(inositol_pat)
    
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    if not inositol_matches:
        return False, "No inositol ring found"
    
    # Create sets of all atom indices that are part of glycerol and inositol matches.
    glycerol_atoms = set()
    for match in glycerol_matches:
        glycerol_atoms.update(match)
    inositol_atoms = set()
    for match in inositol_matches:
        inositol_atoms.update(match)
    
    # Define a SMARTS for a phosphate fragment:
    # [P](=O)(O)(O) matches phosphorus=O and three single-bonded oxygens.
    phosphate_smarts = "[P](=O)(O)(O)"
    phosphate_pat = Chem.MolFromSmarts(phosphate_smarts)
    if phosphate_pat is None:
        return False, "Error in phosphate SMARTS pattern"
    
    phosphate_matches = mol.GetSubstructMatches(phosphate_pat)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Now check for at least one phosphate that has two distinct oxygen substituents:
    # one attached to an atom from the glycerol set and one attached to an atom from the inositol set.
    for match in phosphate_matches:
        # In our SMARTS, the first atom should be the phosphorus.
        # (The ordering of the remaining atoms is not guaranteed so look at the neighbors.)
        # Find the phosphorus atom in this matched fragment by checking the atomic number.
        p_idx = None
        for idx in match:
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 15:
                p_idx = idx
                break
        if p_idx is None:
            continue  # should not happen
        
        P_atom = mol.GetAtomWithIdx(p_idx)
        # Flags to signal if we see the required bonds
        glycerol_link = False
        inositol_link = False
        
        # Loop over all neighbors of the phosphorus.
        for neighbor in P_atom.GetNeighbors():
            # We are interested in oxygen atoms.
            if neighbor.GetAtomicNum() != 8:
                continue
            oxy_idx = neighbor.GetIdx()
            # For this oxygen, check its neighbors (other than the phosphorus).
            for nb in neighbor.GetNeighbors():
                if nb.GetIdx() == p_idx:
                    continue
                # If the neighbor (usually a carbon) belongs to the glycerol fragment match:
                if nb.GetIdx() in glycerol_atoms:
                    glycerol_link = True
                # If it belongs to the inositol fragment match:
                if nb.GetIdx() in inositol_atoms:
                    inositol_link = True
            # If both links are found from this phosphate, we can classify as glycerophosphoinositol.
            if glycerol_link and inositol_link:
                return True, ("Molecule contains a glycerol backbone with a phosphate group that bridges "
                              "to an inositol headgroup")
    
    return False, "Phosphate group is not properly bridging a glycerol and an inositol moiety"

# Example usage:
if __name__ == "__main__":
    # This is a test SMILES string (a simplified example); in practice you would test with a variety of structures.
    test_smiles = "C(COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O)OC(=O)CCCCCCCCCCCCCCC"
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)