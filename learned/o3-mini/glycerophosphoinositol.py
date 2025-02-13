"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: glycerophosphoinositol
Definition: Any glycerophospholipid having the polar alcohol inositol 
esterified to the phosphate group at the sn-3 position of the glycerol backbone.
"""

from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    
    Our strategy is as follows:
      1. Identify a glycerol backbone carrying a phosphate at the sn-3 position.
         We use a specific SMARTS: the first carbon is a non‐ring CH2 with OH,
         the middle carbon is CH with OH, and the third (sn‑3) carbon is a non‐ring CH2
         with a phosphate group attached (–OP(=O)(O)O).
      2. Identify the myo‐inositol headgroup using an explicit cyclic SMARTS.
      3. From the glycerol match, identify the phosphate atom (via the sn‑3 carbon).
         Then, verify that at least one oxygen bound to phosphorus (other than the oxygen 
         linking back to glycerol) connects to a carbon that is a member of an inositol ring.
      4. If found, classification is positive.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a glycerophosphoinositol, else False.
        str: Explanation for the decision.
    """
    
    # Parse the input SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS for myo-inositol:
    # This pattern represents a six-membered cyclohexane ring carrying one hydroxyl on each carbon.
    inositol_smarts = "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    inositol_pat = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pat is None:
        return False, "Error in inositol SMARTS pattern"
    
    # Define a SMARTS for a glycerol backbone with phosphate at sn-3.
    # We require:
    #   - A non-ring primary carbon [CH2;!R] with –OH (sn-1)
    #   - A middle CH (with –OH) (sn-2)
    #   - A non-ring primary carbon [CH2;!R] substituted with a phosphate group (sn-3)
    # The phosphate part is specified as OP(=O)(O)O.
    glycerol_smarts = "[CH2;!R](O)[CH](O)[CH2;!R](OP(=O)(O)O)"
    glycerol_pat = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol_pat is None:
        return False, "Error in glycerol SMARTS pattern"
        
    glycerol_matches = mol.GetSubstructMatches(glycerol_pat)
    if not glycerol_matches:
        return False, "No glycerol-phosphate backbone found"

    # Get inositol matches and build a set of atom indices that are part of any inositol ring.
    inositol_matches = mol.GetSubstructMatches(inositol_pat)
    if not inositol_matches:
        return False, "No inositol ring found"
    inositol_atom_set = set()
    for match in inositol_matches:
        inositol_atom_set.update(match)
    
    # For each glycerol-phosphate match, verify proper connectivity to inositol.
    # Our glycerol pattern is defined with three carbon atoms in order.
    # The third carbon (sn-3) should have a phosphate group.
    for match in glycerol_matches:
        # It is expected that match[0], match[1], match[2] correspond to the 3 glycerol carbons.
        sn3_carbon_idx = match[2]
        sn3_carbon = mol.GetAtomWithIdx(sn3_carbon_idx)
        phosphate_atom = None
        # Identify the phosphorus atom attached to the sn-3 carbon.
        for neighbor in sn3_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 15:  # phosphorus
                phosphate_atom = neighbor
                break
        if phosphate_atom is None:
            continue  # try next match
        
        # Check that phosphate_atom has the expected O=P and three -O substituents.
        # (We do a minimal check looking for at least one double-bonded oxygen.)
        p_neighbors = list(phosphate_atom.GetNeighbors())
        has_doubly_bound_O = any(nb.GetAtomicNum()==8 and 
                                 any(bond.GetBondTypeAsDouble()==True for bond in phosphate_atom.GetBonds() if bond.GetOtherAtom(phosphate_atom)==nb)
                                 for nb in p_neighbors)
        if not has_doubly_bound_O:
            continue  # does not look like a phosphate
        
        # Now, of the oxygen atoms attached to P, one should be already linked to the glycerol sn-3 carbon.
        # Look at the remaining oxygen atoms and see if at least one is connected (via a carbon) to the inositol ring.
        glycerol_o_idx = None
        # Determine which oxygen on phosphate comes from glycerol.
        for nb in phosphate_atom.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                for sub_nb in nb.GetNeighbors():
                    if sub_nb.GetIdx() == sn3_carbon_idx:
                        glycerol_o_idx = nb.GetIdx()
                        break
                if glycerol_o_idx is not None:
                    break

        # Now check the other oxygens bonded to P:
        inositol_link_found = False
        for nb in phosphate_atom.GetNeighbors():
            if nb.GetAtomicNum() != 8:
                continue
            # Skip the oxygen that links back to glycerol.
            if nb.GetIdx() == glycerol_o_idx:
                continue
            # For each such oxygen, check its neighbors (other than phosphate).
            for oxy_nei in nb.GetNeighbors():
                if oxy_nei.GetIdx() == phosphate_atom.GetIdx():
                    continue
                # If any neighbor carbon is part of an inositol match, consider it a link.
                if oxy_nei.GetAtomicNum() == 6 and oxy_nei.GetIdx() in inositol_atom_set:
                    inositol_link_found = True
                    break
            if inositol_link_found:
                break
        
        if inositol_link_found:
            return True, ("Molecule contains a glycerol backbone with a phosphate group at the sn-3 position "
                          "bridging to a myo-inositol headgroup")
    
    return False, "No phosphate group found that directly bridges a defined glycerol and an inositol fragment"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES for glycerophosphoinositol.
    test_smiles = "[C@@H]1([C@@H]([C@@H]([C@@H]([C@H]([C@@H]1O)O)O)O)OP(OC[C@@H](COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O)(=O)O)O"
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)