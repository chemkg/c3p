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
  1. Loop over each phosphorus atom in the molecule.
  2. For that phosphorus (P), check its oxygen neighbors.
      a. At least one oxygen (O_inos) should be connected on its “other side”
         to a carbon that is part of a 6-membered ring (heuristically assigned 
         to an inositol ring).
      b. A different oxygen (O_amide) should be connected to a carbon that is part 
         of an amide group (i.e. the carbon has a double-bonded oxygen and is bonded 
         to a nitrogen).
  3. Also verify that the molecule shows evidence of having long fatty acyl chains 
     (by counting carbon atoms and rotatable bonds).
     
If these criteria are met, we classify the compound as an inositol phosphoceramide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.

    The detection steps are:
      1. Identify phosphorus atoms as candidate bridge centers.
      2. For each phosphorus, require one oxygen neighbor that is linked (via its other
         bond) to a ring member of a 6-membered ring (inositol candidate) and another oxygen
         neighbor that leads to a carbon that is part of an amide (i.e. has a double-bonded oxygen
         and a nitrogen neighbor).
      3. Require that the molecule supports long fatty acyl chains by having a large carbon count
         and many rotatable bonds.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the inositol phosphoceramide criteria, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: check if a carbon is part of an amide group.
    def is_amide_carbon(carbon, mol):
        # For a candidate carbon atom, we check:
        #  1) It must have at least one double-bonded oxygen neighbor.
        #  2) It must have a neighboring nitrogen (amide N).
        has_dbl_o = False
        has_nitrogen = False
        for nbr in carbon.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                has_dbl_o = True
            if nbr.GetAtomicNum() == 7:
                has_nitrogen = True
        return has_dbl_o and has_nitrogen

    bridge_found = False
    # Loop over candidates of phosphorus atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue  # not phosphorus
        # Get oxygen neighbors of phosphorus
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(oxy_neighbors) < 2:
            continue  # not enough oxygen bridges
        
        found_inositol = False
        found_amide = False
        # Check for the inositol side:
        # We require that at least one oxygen neighbor of P is bonded (besides P) 
        # to an atom that is in a 6-membered ring. We assume the inositol ring is 6-membered.
        for oxy in oxy_neighbors:
            for nbr in oxy.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                # Check if the neighbor atom is part of a 6-membered ring.
                if nbr.IsInRingSize(6):
                    found_inositol = True
                    break
            if found_inositol:
                break
        
        # Now, check for the ceramide/amidic side:
        # Look for an oxygen (other than the one we used for the inositol connection)
        # that leads to a carbon that is part of an amide (C(=O)N)
        for oxy in oxy_neighbors:
            for nbr in oxy.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:  # potential carbon linking to the ceramide chain
                    if is_amide_carbon(nbr, mol):
                        found_amide = True
                        break
            if found_amide:
                break

        if found_inositol and found_amide:
            bridge_found = True
            break

    if not bridge_found:
        return False, "Phosphodiester bridge linking inositol ring and ceramide moiety not detected"

    # Optional step: Check for long aliphatic chains.
    # We ensure that there are a sufficient number of carbons and many rotatable bonds.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 20:
        return False, "Too few carbons to support long acyl chains in a ceramide moiety"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds; fatty acyl chains may be too short"
    
    return True, "Molecule contains a bridging phosphate with an inositol-linked oxygen and an amide-linked carbon (indicative of a ceramide), plus long aliphatic chains."

# Example usage:
if __name__ == "__main__":
    # Example SMILES (one of the provided compounds)
    test_smiles = ("CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](COP(O)(=O)O"
                   "[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)"
                   "[C@H](O)C(O)CCCCCCCCCCCCCC")
    result, reason = is_inositol_phosphoceramide(test_smiles)
    print(result, reason)