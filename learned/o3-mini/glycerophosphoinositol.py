"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: glycerophosphoinositol
Definition: Any glycerophospholipid having the polar alcohol inositol 
esterified to the phosphate group at the sn-3 position of the glycerol backbone.

Improvement notes:
 - Use a more specific SMARTS pattern for glycerol (3‐carbon chain with two primary alcohols)
 - Use an inositol SMARTS for a myo‐inositol ring
 - When examining phosphate groups ([P](=O)(O)(O)), require that one oxygen substituent
   is connected to a carbon that is uniquely part of the glycerol fragment and another oxygen
   is connected to a carbon that is uniquely part of the inositol ring.
"""

from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    
    The strategy is to (1) identify a simplified glycerol backbone and an inositol ring;
    (2) find a phosphate ([P](=O)(O)(O)) group; and (3) check if that phosphate has two
    oxygen substituents such that one is linked (through a carbon) to a unique glycerol atom
    and a second is linked to a unique inositol atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a glycerophosphoinositol, else False.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns.
    # A simplified myo-inositol pattern: six-membered ring with one hydroxyl on each carbon.
    # (This will match a cyclohexane in which every carbon is substituted with at least one OH.)
    inositol_smarts = "C1C(O)C(O)C(O)C(O)C1O"
    inositol_pat = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pat is None:
        return False, "Error in inositol SMARTS pattern"
    
    # Define a simplified glycerol backbone pattern.
    # Here we want a three-carbon chain with terminal primary hydroxyl groups.
    # One common representation is: HO-C(CO)CO which covers that the terminal carbons 
    # (sn-1 and sn-3) are CH2OH and the middle one is CHOH.
    glycerol_smarts = "[#6;X4](O)[#6;X4](O)[#6;X4](O)"
    glycerol_pat = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol_pat is None:
        return False, "Error in glycerol SMARTS pattern"
    
    # Find all matches for inositol and glycerol fragments.
    inositol_matches = mol.GetSubstructMatches(inositol_pat)
    glycerol_matches = mol.GetSubstructMatches(glycerol_pat)
    
    if not inositol_matches:
        return False, "No inositol ring found"
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    
    # Build sets of atom indices for each fragment.
    glycerol_atoms = set()
    for match in glycerol_matches:
        glycerol_atoms.update(match)
    inositol_atoms = set()
    for match in inositol_matches:
        inositol_atoms.update(match)
    
    # Remove atoms that occur in both sets (to avoid overlap that causes misclassification).
    unique_glycerol_atoms = glycerol_atoms.difference(inositol_atoms)
    unique_inositol_atoms = inositol_atoms.difference(glycerol_atoms)
    if not unique_glycerol_atoms:
        return False, "Glycerol fragment not clearly defined (overlap with inositol?)"
    if not unique_inositol_atoms:
        return False, "Inositol fragment not clearly defined (overlap with glycerol?)"
    
    # Define phosphate SMARTS: phosphorus with one double bond O and three single-bonded oxygens.
    phosphate_smarts = "[P](=O)(O)(O)"
    phosphate_pat = Chem.MolFromSmarts(phosphate_smarts)
    if phosphate_pat is None:
        return False, "Error in phosphate SMARTS pattern"
    
    phosphate_matches = mol.GetSubstructMatches(phosphate_pat)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # For each phosphate group found, check its neighbors:
    for match in phosphate_matches:
        p_atom = None
        # Identify the phosphorus atom within the match (atomic number 15)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 15:
                p_atom = atom
                break
        if p_atom is None:
            continue  # safety check
        
        glycerol_link_found = False
        inositol_link_found = False
        
        # Look at neighbors of phosphorus
        for neighbor in p_atom.GetNeighbors():
            # We want oxygen neighbors (atomic number 8)
            if neighbor.GetAtomicNum() != 8:
                continue
            oxy_idx = neighbor.GetIdx()
            # For each oxygen neighbor, check its other neighbors (should be a carbon)
            for nb in neighbor.GetNeighbors():
                # Skip the phosphorus atom itself.
                if nb.GetIdx() == p_atom.GetIdx():
                    continue
                # We require that the neighbor is a carbon
                if nb.GetAtomicNum() != 6:
                    continue
                c_idx = nb.GetIdx()
                # If this carbon is uniquely in the glycerol backbone, mark glycerol linkage.
                if c_idx in unique_glycerol_atoms:
                    glycerol_link_found = True
                # If this carbon is uniquely in the inositol ring, mark inositol linkage.
                if c_idx in unique_inositol_atoms:
                    inositol_link_found = True
            # If both links are already found, we can classify the phosphate as bridging.
            if glycerol_link_found and inositol_link_found:
                return True, ("Molecule contains a glycerol backbone with a phosphate group that bridges "
                              "to an inositol headgroup")
    
    return False, "No phosphate group found that directly bridges a unique glycerol and a unique inositol fragment"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES for a glycerophosphoinositol.
    test_smiles = "[C@@H]1([C@@H]([C@@H]([C@@H]([C@H]([C@@H]1O)O)O)O)OP(OC[C@@H](COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O)(=O)O)O"
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)