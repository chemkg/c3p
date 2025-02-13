"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: Carbonate Ester

Definition:
  A carbonate ester is any derivative of carbonic acid (O=C(OH)2) in which one or both –OH groups
  have been replaced by organyl groups (R). This includes diesters (O=C(–O–R)2), monoesters (O=C(–OH)(–O–R))
  and cyclic carbonates that embed the C(=O)(O)(O) unit.
  
The strategy is to:
  1. Parse the SMILES and add explicit hydrogens so that –OH groups will be represented as O–H.
  2. Loop over all carbon atoms and identify those that have exactly three heavy (non‐hydrogen) neighbors,
     all being oxygen.
  3. Among its bonds, require exactly one double bond (i.e. a carbonyl) and two single bonds to oxygen.
  4. For each of the single‐bonded oxygens, examine if it is substituted (i.e. if it has any heavy neighbor 
     in addition to the central carbon). In a valid carbonate ester, at least one oxygen must be substituted.
  
If such a carbonate moiety is found, the function returns True along with a reason indicating whether
one or both –OH groups have been substituted; otherwise it returns False.
"""

from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    
    The method adds explicit hydrogens so that –OH groups are represented properly.
    Then for every carbon atom in the molecule, it checks whether the carbon appears to be in a
    carbonic acid environment (C with one double-bonded oxygen and two single-bonded oxygens), and further
    whether at least one of the two single-bonded oxygens is substituted (i.e. attached to a non-hydrogen besides the carbon).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a carbonate ester, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups are seen as O-H.
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms to try to identify the carbonate moiety
    for atom in mol.GetAtoms():
        # We only focus on carbon atoms
        if atom.GetAtomicNum() != 6:
            continue
        
        # Get non-hydrogen (heavy) neighbors
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # The central carbon in a carbonic acid unit should be bound to exactly 3 heavy atoms (all oxygens)
        if len(heavy_neighbors) != 3:
            continue
        
        # Check that all three neighbors are oxygens.
        if any(nbr.GetAtomicNum() != 8 for nbr in heavy_neighbors):
            continue
        
        # Now check the bond types: exactly one double bond (C=O) and two single bonds (C-O)
        double_bonded_oxygen = None
        single_bonded_oxygens = []
        for bond in atom.GetBonds():
            # Consider only bonds where the neighbor is oxygen
            # Note: We check both ends but we want the oxygen neighbor.
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() != 8:
                continue
            # Distinguish double and single bonds
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                if double_bonded_oxygen is None:
                    double_bonded_oxygen = nbr
                else:
                    # More than one double bond disqualifies this carbon.
                    double_bonded_oxygen = None
                    break
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                single_bonded_oxygens.append(nbr)
        # Must have exactly one double-bonded oxygen and exactly two single-bonded oxygens.
        if double_bonded_oxygen is None or len(single_bonded_oxygens) != 2:
            continue
        
        # Check substitution on the two single-bonded oxygens:
        # For each such oxygen, consider its heavy neighbors excluding the central carbon.
        substitution_flags = []
        for oxy in single_bonded_oxygens:
            other_heavy = [nbr for nbr in oxy.GetNeighbors() if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() > 1]
            # In an unsubstituted –OH group, the only heavy neighbor would be the central carbon.
            # So if other_heavy is non-empty, the O has an extra (organyl) substituent.
            substitution_flags.append(len(other_heavy) > 0)
        
        # At least one oxygen must be substituted (for a carbonate ester).
        if any(substitution_flags):
            if all(substitution_flags):
                reason = "Found carbonate moiety (O=C(–O–R)2) with both –OH replaced by organyl groups."
            else:
                reason = "Found carbonate moiety (O=C(–OH)(–O–R)) with at least one organyl substitution."
            return True, reason

    # If no such substructure is found, then we do not classify the molecule as a carbonate ester.
    return False, "No carbonate ester motif (C(=O)(O)(O) derivative with required substitution) detected."

# Uncomment below for quick testing:
# test_smiles = [
#     "O(CCC(C)=C)C(OCC)=O",  # Ethyl 3-methylbut-3-enyl carbonate (diester)
#     "COC(=O)OC",            # dimethyl carbonate (diester)
#     "O1CCOC1=O",            # Ethylene carbonate (cyclic diester)
#     "CCOC(O)=O",            # monoethyl carbonate (monoester)
#     "CC(C)(C)C(=O)OCOP(=O)(COCCN1C=NC2=C1N=CN=C2N)OCOC(=O)C(C)(C)C",  # adefovir pivoxil (previously a false negative)
#     "O=C1OC(=O)C=CO1",       # Pleurone (cyclic carbonate, previously a false negative)
#     "CCOc1nc2cccc(C(=O)OCc3oc(=O)oc3C)c2n1Cc1ccc(cc1)-c1ccccc1-c1noc(=O)[nH]1",  # azilsartan medoxomil (previously a false negative)
# ]
# for smi in test_smiles:
#     result, reason = is_carbonate_ester(smi)
#     print(f"SMILES: {smi}\n --> {result}. Reason: {reason}\n")