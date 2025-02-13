"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid
Definition: A derivative of glycerol in which one hydroxy group (commonly but not necessarily primary)
is esterified with phosphoric acid and the other two are esterified with fatty acids.
Our revised approach (approximate):
  1. The molecule must be valid.
  2. It must contain exactly one phosphorus atom.
  3. It should not contain any nitrogen atoms.
  4. It must not contain any rings (to favor an acyclic glycerol backbone).
  5. The phosphate group should have exactly 3 oxygen neighbors (helping to exclude entities with an extra head‐group).
  6. The molecule must show exactly two fatty acid ester groups. We define a fatty acyl ester as an oxygen
     (not directly bonded to phosphorus) connected to a carbonyl – i.e. the pattern “[O;!$([O]-[#15])]-C(=O)[#6]”.
  7. Finally, the phosphate ester linkage is checked by confirming that at least one oxygen attached to the phosphorus
     is connected to an sp3 carbon that is not part of a carbonyl.
If all criteria are met we return True.
Note: This “rule‐based” classifier is approximate.
"""

from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.

    The criteria (approximate) are:
      - Valid molecule.
      - Exactly one phosphorus (P) atom.
      - No nitrogen atoms are allowed.
      - No rings (the glycerol backbone is expected to be acyclic).
      - The phosphate group must have exactly 3 oxygen neighbors,
        meaning it is not extended into an extra head‐group (e.g. glycerol as in PG).
      - Exactly 2 fatty acid ester groups are detected via the SMARTS
        pattern "[O;!$([O]-[#15])]-C(=O)[#6]".
      - At least one phosphate ester linkage is found (an oxygen on P is bonded to an sp3 carbon that is not carbonylated).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a phosphatidic acid, else False.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Criterion 1: The molecule must contain exactly one phosphorus atom.
    phos_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phos_atoms) != 1:
        return False, f"Expected exactly 1 phosphorus atom, found {len(phos_atoms)}"
    p_atom = phos_atoms[0]

    # Criterion 2: Molecule must not contain nitrogen.
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Presence of nitrogen detected – likely not phosphatidic acid"

    # Criterion 3: Molecule should not contain any rings. (Helps to exclude extra head‐groups.)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings – expected an acyclic glycerol backbone"

    # Criterion 4: For PA the phosphate group should be simple.
    # Check that phosphorus has exactly three oxygen neighbors.
    p_oxy_neighbors = [nbr for nbr in p_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(p_oxy_neighbors) != 3:
        return False, f"Expected phosphate to have exactly 3 oxygen neighbors, found {len(p_oxy_neighbors)}"

    # Criterion 5: Look for fatty acid ester groups.
    # Our SMARTS matches an oxygen not directly bound to phosphorus, directly linked to a carbonyl: -O-C(=O)-R.
    fatty_ester_smarts = "[O;!$([O]-[#15])]-C(=O)[#6]"
    fatty_ester_pattern = Chem.MolFromSmarts(fatty_ester_smarts)
    matches = mol.GetSubstructMatches(fatty_ester_pattern)
    fatty_ester_oxygens = set(match[0] for match in matches)
    if len(fatty_ester_oxygens) != 2:
        return False, f"Expected exactly 2 fatty acid ester groups, found {len(fatty_ester_oxygens)}"

    # Criterion 6: Check for a phosphate ester linkage.
    # We require at least one oxygen attached to phosphorus to be linked to an sp3 (non-carbonyl) carbon.
    phosphate_linkage_found = False
    for o_atom in p_oxy_neighbors:
        # Look at neighbors of the oxygen (ignoring the phosphorus)
        for subnbr in o_atom.GetNeighbors():
            if subnbr.GetIdx() == p_atom.GetIdx():
                continue
            if subnbr.GetAtomicNum() == 6:  # a carbon atom
                # Check if the carbon is not carbonylated,
                # i.e. it does not have a double bond to an oxygen.
                is_carbonyl = False
                for bond in subnbr.GetBonds():
                    # Using bond type as a proxy for double bond
                    if bond.GetBondType().name == "DOUBLE":
                        other = bond.GetOtherAtom(subnbr)
                        if other.GetAtomicNum() == 8:
                            is_carbonyl = True
                            break
                if not is_carbonyl:
                    phosphate_linkage_found = True
                    break
        if phosphate_linkage_found:
            break
    if not phosphate_linkage_found:
        return False, "No phosphate ester linkage to a non-carbonyl (glycerol) carbon found"

    return True, ("Molecule has one phosphorus, no nitrogen, an acyclic backbone, exactly 3 oxygens attached to the phosphorus, "
                  "exactly 2 fatty acid ester groups, and a proper phosphate ester linkage – consistent with phosphatidic acid.")

# For testing purposes, one might try calling the function with one of the provided SMILES.
if __name__ == "__main__":
    # Example: PA(14:0/16:1(9Z))
    test_smiles = "[C@](COC(=O)CCCCCCCCCCCCC)(OC(=O)CCCCCCC/C=C\\CCCCCC)([H])COP(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)