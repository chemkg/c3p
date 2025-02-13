"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid
Definition: A derivative of glycerol in which one hydroxy group (commonly but not necessarily primary) 
is esterified with phosphoric acid and the other two are esterified with fatty acids.
Our approach (approximate):
  1. The molecule must be valid.
  2. It must contain exactly one phosphorus atom.
  3. It should not contain nitrogen or rings (to help exclude other phospholipids with extra head‐groups).
  4. It must show at least two fatty acid ester groups – defined by the SMARTS pattern [O]-C(=O)[#6] 
     (but we ignore any O directly attached to phosphorus).
  5. It must show a phosphate ester linkage – we check that at least one oxygen attached to the unique 
     phosphorus is itself bound to an sp3 carbon that does not appear to be carbonylated.
If these criteria are met we return True.
Note: This “rule‐based” classifier is approximate.
"""

from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    
    Our criteria (approximate):
      - Valid molecule.
      - Exactly one phosphorus (P) atom.
      - No nitrogen atoms (to exclude headgroups from e.g. PC, PE, etc).
      - No rings (PA backbone should be acyclic).
      - At least two fatty acid ester groups (–O–C(=O)–R) not directly attached to P.
      - Presence of a phosphate ester linkage: at least one oxygen bound to the P
          is also bound to an aliphatic carbon (i.e. not a carbonyl carbon).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Exactly one phosphorus atom (atomic number 15)
    phos_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phos_atoms) != 1:
        return False, f"Expected exactly 1 phosphorus atom, found {len(phos_atoms)}"
        
    # Criterion 2: Exclude molecules with nitrogen (N) since extra headgroups (like choline, ethanolamine, etc) are present in other lipids.
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Presence of nitrogen detected – likely not phosphatidic acid"
    
    # Criterion 3: Exclude molecules containing rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings – expected an acyclic glycerol backbone"
    
    # Criterion 4: Look for fatty acid ester groups.
    # A fatty acid ester group (as attached from the glycerol backbone) is approximated by the pattern:
    #   [O]-C(=O)[#6]
    # We ignore any match where the oxygen is directly attached to a phosphorus.
    fatty_ester_pattern = Chem.MolFromSmarts("[O]-C(=O)[#6]")
    matches = mol.GetSubstructMatches(fatty_ester_pattern)
    fatty_ester_oxygens = set()
    for match in matches:
        # 'match[0]' corresponds to the oxygen atom in the pattern.
        oxy_atom = mol.GetAtomWithIdx(match[0])
        # Exclude if any neighbor of this oxygen is phosphorus.
        if any(neighbor.GetAtomicNum() == 15 for neighbor in oxy_atom.GetNeighbors()):
            continue
        fatty_ester_oxygens.add(match[0])
    if len(fatty_ester_oxygens) < 2:
        return False, f"Expected at least 2 fatty acid ester groups, found {len(fatty_ester_oxygens)}"
    
    # Criterion 5: Check for a phosphate ester linkage.
    # We require that at least one oxygen attached to the unique phosphorus is also linked to an sp3 carbon
    # that does not appear to be a carbonyl carbon.
    p_atom = phos_atoms[0]
    phosphate_linkage_found = False
    for nbr in p_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 8:  # oxygen neighbor
            # Look at each neighbor of this oxygen (other than the phosphorus)
            for subnbr in nbr.GetNeighbors():
                if subnbr.GetIdx() == p_atom.GetIdx():
                    continue
                # We expect a carbon (atomic num 6) that is sp3 rather than double-bonded to oxygen (carbonyl).
                if subnbr.GetAtomicNum() == 6:
                    # Check bonds of this carbon: if any bond (other than to oxygen nbr) is a double bond to O, skip.
                    is_carbonyl = False
                    for bond in subnbr.GetBonds():
                        if bond.GetBondType().name == "DOUBLE":
                            # If the other atom is oxygen then we consider this carbonyl (thus not the glycerol linkage).
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
        return False, "No phosphate ester linkage found (an oxygen on P is not linked to an aliphatic carbon)"
    
    # If all tests pass, assume the molecule is phosphatidic acid.
    return True, ("Molecule has one phosphate, at least two fatty acid ester groups and a phosphate "
                   "ester linkage with no extra headgroup features (e.g. nitrogen or rings), consistent with phosphatidic acid.")

# For testing purposes (one can call the function with one of the provided SMILES):
if __name__ == "__main__":
    # Example: PA(14:0/16:1(9Z))
    test_smiles = "[C@](COC(=O)CCCCCCCCCCCCC)(OC(=O)CCCCCCC/C=C\\CCCCCC)([H])COP(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)