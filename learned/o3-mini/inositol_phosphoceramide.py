"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
#!/usr/bin/env python
"""
Classifies: inositol phosphoceramide
Definition:
  A phosphosphingolipid in which an inositol residue and the ceramide moiety 
  are linked via a phosphodiester bridge. The ceramide moiety contains 
  substituents R(1) and R(2) which vary with different sphingoid bases and 
  fatty acyl moieties.
  
Detection strategy:
  1. Use SMARTS queries to detect an inositol substructure (a cyclohexanol with 6 hydroxyl groups)
     and an amide function (as indicative of the ceramide part).
  2. Identify phosphorus atoms and for each such atom, check:
       a. At least one oxygen neighbor is part of an inositol match.
       b. At least one oxygen neighbor (via its bonded atom) is part of an amide match.
  3. Verify that the molecule has many carbons and a fair number of rotatable bonds 
     (to ensure the presence of long aliphatic chains).
  
If these criteria are met, the compound is classified as an inositol phosphoceramide.
 
Note: The SMARTS for the inositol pattern is specific and may not catch every variant,
but it seems to match many of the examples provided.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    Detection steps:
      1. Pre-match an inositol group using a SMARTS pattern typical for cyclohexane-1,2,3,4,5,6-hexol.
      2. Pre-match an amide motif (C(=O)N) considered part of the ceramide moiety.
      3. Loop over every phosphorus atom in the molecule. For each phosphorus, look at its oxygen neighbors:
             a. One of the oxygens must belong to an inositol substructure.
             b. One oxygen (or another neighboring atom in the oxygen’s bonds) must be part of an amide motif.
      4. Also check that there are enough carbon atoms and rotatable bonds to support long acyl chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the inositol phosphoceramide criteria, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Pre-match inositol: a rough SMARTS pattern for an inositol ring (cyclohexanol with 6 OH groups).
    # This pattern covers a typical arrangement: O-[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O
    inositol_smarts = "O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"
    inositol_query = Chem.MolFromSmarts(inositol_smarts)
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "No inositol substructure detected"
    
    # Pre-match the amide motif, indicative of the ceramide moiety.
    amide_smarts = "C(=O)N"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    amide_matches = mol.GetSubstructMatches(amide_query)
    if not amide_matches:
        return False, "No amide substructure (ceramide indication) detected"
    
    # Optional: Ensure the molecule supports long fatty acyl chains.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 20:
        return False, "Too few carbons to support long acyl chains in the ceramide moiety"
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds; fatty acyl chains may be too short"
    
    # Now, loop over phosphorus atoms to check for the phosphodiester bridge.
    bridge_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:  # Must be phosphorus
            continue
        # Get oxygen neighbors of phosphorus.
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(oxy_neighbors) < 2:
            continue  # Not enough oxygen bridges
      
        inositol_connected = False
        amide_connected = False
        
        # Check that at least one oxygen is directly part of an inositol substructure.
        for oxy in oxy_neighbors:
            for match in inositol_matches:
                if oxy.GetIdx() in match:
                    inositol_connected = True
                    break
            if inositol_connected:
                break

        # Check that at least one oxygen neighbor (via one of its attached atoms) is connected to an amide.
        for oxy in oxy_neighbors:
            for nbr in oxy.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue  # skip back to phosphorus
                for match in amide_matches:
                    if nbr.GetIdx() in match:
                        amide_connected = True
                        break
                if amide_connected:
                    break
            if amide_connected:
                break
        
        if inositol_connected and amide_connected:
            bridge_found = True
            break

    if not bridge_found:
        return False, "Phosphodiester bridge linking inositol ring and ceramide moiety not detected"
    
    return True, "Molecule contains an inositol substructure, an amide-bearing ceramide moiety, and a bridging phosphate."

# Example usage
if __name__ == "__main__":
    # One of the provided example SMILES strings (Man-beta1-6-Ins-1-P-Cer(t18:0/2-OH-24:0))
    test_smiles = ("CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](COP(O)(=O)O"
                   "[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)"
                   "[C@H](O)C(O)CCCCCCCCCCCCCC")
    result, reason = is_inositol_phosphoceramide(test_smiles)
    print(result, reason)