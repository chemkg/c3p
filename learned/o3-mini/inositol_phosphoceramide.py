"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: inositol phosphoceramide
Definition:
  A phosphosphingolipid in which an inositol residue and the ceramide moiety are linked via a phosphodiester bridge.
  The ceramide moiety contains substituents R(1) and R(2) which vary with different sphingoid bases and fatty acyl moieties.

Our improved approach:
  1. Look for an inositol ring (using a relatively generic ring pattern for myo‐inositol,
     ignoring stereochemistry).
  2. For at least one oxygen in the inositol ring, check if it is bonded to a phosphorus atom.
  3. From that phosphorus, follow a non‐inositol oxygen that bonds to a carbon
     which is part of an amide bond (i.e. has a adjacent double‐bonded oxygen and a neighbor nitrogen).
  4. Finally, require that the molecule shows evidence of long fatty acyl chains
     (sufficient carbon count and many rotatable bonds).

If these criteria are met, we classify the compound as an inositol phosphoceramide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.

    The detection steps are:
      1. Detect an inositol ring (neutral, six-membered ring featuring five hydroxyl groups and one ring oxygen).
      2. Check that at least one oxygen on that ring is bonded to a phosphorus atom.
      3. Verify that from that phosphorus a different oxygen leads to a carbon which is part of a C(=O)N (amide) group.
      4. (Optional) Check overall carbon count and number of rotatable bonds to ensure the presence of long fatty acyl chains.
     
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an inositol phosphoceramide, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Look for inositol ring.
    # SMARTS for a generic inositol ring (ignoring stereochemistry):
    inositol_smarts = "OC1C(O)C(O)C(O)C(O)C1O"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring detected"
    
    # Step 2 & 3: For atoms in the inositol ring, look for an oxygen connected to phosphorus,
    # and from that phosphorus look for an oxygen that connects to a carbon involved in an amide.
    bridge_found = False
    # Loop over each inositol match (each is a tuple of atom indices from the pattern)
    for match in inositol_matches:
        # Get the atoms corresponding to the inositol ring.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # We are only interested in oxygen atoms.
            if atom.GetAtomicNum() != 8:
                continue
            # Look at neighbors of this oxygen
            for nbr in atom.GetNeighbors():
                # Check if neighbor is phosphorus (atomic number 15)
                if nbr.GetAtomicNum() == 15:
                    P_atom = nbr
                    # Get oxygen neighbors of the phosphorus
                    p_oxy_neighbors = [a for a in P_atom.GetNeighbors() if a.GetAtomicNum() == 8]
                    # There should be at least two oxygens (one is from inositol, one connecting to the other side)
                    if len(p_oxy_neighbors) < 2:
                        continue
                    # Now, check those oxygens except the one from the inositol ring.
                    for oxy in p_oxy_neighbors:
                        if oxy.GetIdx() == atom.GetIdx():
                            continue
                        # Look for a carbon neighbor of this oxygen (this is our candidate linker)
                        for c in oxy.GetNeighbors():
                            if c.GetAtomicNum() != 6:
                                continue
                            # Check if the carbon is part of an amide group
                            # (i.e., it has a neighbor oxygen with a double bond and a neighboring nitrogen)
                            amide_flag = False
                            for neighbor in c.GetNeighbors():
                                # Look for oxygen double-bonded to this carbon
                                if neighbor.GetAtomicNum() == 8:
                                    bond = mol.GetBondBetweenAtoms(c.GetIdx(), neighbor.GetIdx())
                                    if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                                        # Now check if the carbon also has a neighbor nitrogen (amide N)
                                        for n in c.GetNeighbors():
                                            if n.GetAtomicNum() == 7:
                                                amide_flag = True
                                                break
                                if amide_flag:
                                    break
                            if amide_flag:
                                bridge_found = True
                                break
                        if bridge_found:
                            break
                    if bridge_found:
                        break
            if bridge_found:
                break
        if bridge_found:
            break

    if not bridge_found:
        return False, "Phosphodiester bridge linking inositol ring and ceramide moiety not detected"
    
    # Step 4: Optional checks for long aliphatic chains (fatty acyl chains typically found in ceramides).
    # Check total carbon count.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 20:
        return False, "Too few carbons to support long acyl chains in a ceramide moiety"
    
    # Check for a sufficient number of rotatable bonds.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds; fatty acyl chains may be too short"
    
    return True, "Molecule contains an inositol ring linked (via a phosphodiester bridge) to a ceramide moiety with long acyl chains"

# Example usage:
if __name__ == "__main__":
    # You can test with one of the provided SMILES:
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](COP(O)(=O)O[C@@H]1C(O)C(O)C(O)C(O)C1O)[C@H](O)C(O)CCCCCCCCCCCCCC"
    result, reason = is_inositol_phosphoceramide(test_smiles)
    print(result, reason)