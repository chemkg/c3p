"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: Carbonate Ester
Definition: Any carbonate that is carbonic acid (O=C(OH)2) in which one or both -OH groups 
have been replaced by organyl (non-hydrogen) groups. This includes derivatives that are diesters,
monoesters, or cyclic carbonates.
"""

from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    
    The method looks for a carbon atom that is in the formal environment of carbonic acid:
      - It is double bonded to one oxygen.
      - It is single bonded to two other oxygen atoms.
      - Its only heavy-atom (non-H) neighbors are these three oxygens.
      
    Then for the two single-bonded oxygens we check the substitution:
      - In a parent carbonic acid (O=C(OH)2) both oxygens would be unsubstituted.
      - In a carbonate ester at least one of these oxygens has an additional (non-H) substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a carbonate ester, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (Optional) We could add hydrogens explicitly if needed using Chem.AddHs(mol), 
    # but here we can simply count heavy atoms (atomic number > 1).

    # Loop over all carbon atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # skip non-carbon atoms

        # Get heavy neighbors (atoms with atomic number > 1)
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
        # For a carbonic acid core, the carbon should only be bonded to three heavy atoms.
        if len(heavy_neighbors) != 3:
            continue

        # All heavy neighbors must be oxygen.
        if any(n.GetAtomicNum() != 8 for n in heavy_neighbors):
            continue

        # Now, identify bonds: we need exactly one double bond (C=O) 
        # and two single bonds (C–O) from this carbon.
        double_bonded_oxygen = None
        single_bonded_oxygens = []
        for bond in atom.GetBonds():
            # Get the neighbor atom in this bond:
            nbr = bond.GetOtherAtom(atom)
            # Only consider oxygen neighbors
            if nbr.GetAtomicNum() != 8:
                continue
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                if double_bonded_oxygen is None:
                    double_bonded_oxygen = nbr
                else:
                    # More than one double bond to oxygen? Not our target.
                    double_bonded_oxygen = None
                    break
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                single_bonded_oxygens.append(nbr)
        # Must have exactly one carbonyl oxygen and exactly two singly-bound oxygens.
        if double_bonded_oxygen is None or len(single_bonded_oxygens) != 2:
            continue

        # At this point, atom is in a carbonic acid-like environment: O=C(O)(O)
        # Now check the substitution on the two single-bonded oxygens.
        # In an unsubstituted carbonic acid, each oxygen only has the central C as a heavy neighbor.
        substituted_flags = []
        for oxy in single_bonded_oxygens:
            # Get heavy neighbors of the oxygen other than the central carbon.
            other_heavy = [n for n in oxy.GetNeighbors() if n.GetIdx() != atom.GetIdx() and n.GetAtomicNum() > 1]
            # If there is any other heavy neighbor, then that -OH has been replaced by an organyl group.
            substituted_flags.append(len(other_heavy) > 0)
        
        # For a carbonate ester at least one of the single-bonded oxygens must be substituted.
        if any(substituted_flags):
            # Determine a reason string based on whether one or both oxygens are substituted.
            if all(substituted_flags):
                reason = ("Found carbonate moiety (O=C(–O–R)2) with both –OH replaced by organyl groups.")
            else:
                reason = ("Found carbonate moiety (O=C(–OH)(–O–R)) with at least one organyl substitution.")
            return True, reason

    # No valid carbonate ester substructure was found.
    return False, "No carbonate ester motif (O=C(OH)2 derivative with required substitutions) detected."

# Uncomment below lines for quick testing.
# test_smiles = [
#     "O(CCC(C)=C)C(OCC)=O",  # Ethyl 3-methylbut-3-enyl carbonate
#     "COC(=O)OC",            # dimethyl carbonate
#     "O1CCOC1=O",            # Ethylene carbonate
#     "CCOC(O)=O",            # monoethyl carbonate
#     "O=C1OC(=O)C=CO1",       # Pleurone (cyclic carbonate)
#     "CC(C)(C)C(=O)OCOP(=O)(COCCN1C=NC2=C1N=CN=C2N)OCOC(=O)C(C)(C)C",  # adefovir pivoxil
#     "CCOc1nc2cccc(C(=O)OCc3oc(=O)oc3C)c2n1Cc1ccc(cc1)-c1ccccc1-c1noc(=O)[nH]1",  # azilsartan medoxomil
# ]
# for smi in test_smiles:
#     result, reason = is_carbonate_ester(smi)
#     print(f"SMILES: {smi} \n --> {result}. Reason: {reason}\n")