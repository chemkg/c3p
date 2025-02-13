"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid

Definition: A derivative of glycerol in which one hydroxy group – commonly but not necessarily primary – 
is esterified with phosphoric acid and the other two are esterified with fatty acids.

Revised approach (approximate):
  1. The molecule must be valid.
  2. It must contain exactly one phosphorus atom.
  3. It should not contain any nitrogen atoms.
  4. It must have exactly two fatty acid ester groups.
     (A fatty acid ester is defined as an oxygen (not directly bound to phosphorus)
      connected to a carbonyl group via the pattern "[O;!$([O]-[#15])]-C(=O)[#6]".)
  5. In addition, the molecule must contain a glycerol backbone defined as a contiguous three‐carbon chain 
     in which one of the –OH substituents is a phosphate ester (an oxygen attached to P) 
     and the other two substituents are fatty acid ester linkages (an oxygen attached to a carbonyl).
     
If all criteria are met we return True with a reason; else, False with an error explanation.
"""

from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    
    Steps:
      - Valid molecule.
      - Exactly one phosphorus.
      - No nitrogen.
      - Exactly 2 fatty acid ester groups (using SMARTS).
      - Presence of a glycerol backbone: a contiguous 3-carbon chain where 
        one carbon is substituted by an oxygen attached to phosphorus and 
        the other two are substituted by oxygens that are in fatty acid ester linkages.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is phosphatidic acid, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Criterion 1: exactly one phosphorus atom
    phos_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phos_atoms) != 1:
        return False, f"Expected exactly 1 phosphorus atom, found {len(phos_atoms)}"
        
    # Criterion 2: no nitrogen
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Presence of nitrogen detected – not a phosphatidic acid"
        
    # Criterion 3: Check for exactly 2 fatty acid ester groups.
    # We use the SMARTS: oxygen NOT bound to phosphorus, directly attached to a carbonyl.
    fatty_ester_smarts = "[O;!$([O]-[#15])]-C(=O)[#6]"
    fatty_ester_pattern = Chem.MolFromSmarts(fatty_ester_smarts)
    matches = mol.GetSubstructMatches(fatty_ester_pattern)
    # Collect the oxygen atom indices that are the first atom in the match, to avoid double counting.
    fatty_ester_oxygens = {match[0] for match in matches}
    if len(fatty_ester_oxygens) != 2:
        return False, f"Expected exactly 2 fatty acid ester groups, found {len(fatty_ester_oxygens)}"
        
    # Criterion 4: Check for a glycerol backbone.
    # We require the presence of a contiguous three‐carbon chain (glycerol) in which:
    #   – one of the carbons has an oxygen substituent that, in turn, is linked to phosphorus (phosphate ester),
    #   – and the other two carbons have oxygen substituents that are fatty acid ester oxygens.
    
    def is_phosphate_oxygen(o_atom, exclude_idx=None):
        # Returns True if this oxygen is connected (besides possibly to an atom we want to exclude)
        # to a phosphorus atom.
        for nbr in o_atom.GetNeighbors():
            if exclude_idx is not None and nbr.GetIdx() == exclude_idx:
                continue
            if nbr.GetAtomicNum() == 15:
                return True
        return False

    def is_fatty_ester_oxygen(o_atom, exclude_idx=None):
        # Returns True if this oxygen is attached (besides possibly to an atom we want to exclude)
        # to a carbon that is double-bonded to an oxygen (C=O).
        for nbr in o_atom.GetNeighbors():
            if exclude_idx is not None and nbr.GetIdx() == exclude_idx:
                continue
            if nbr.GetAtomicNum() == 6:  # carbon
                # Check if this carbon nbr bears a double bond to an oxygen (a carbonyl)
                for bond in nbr.GetBonds():
                    # If bond is double and the other atom is oxygen, count as part of a C=O
                    if bond.GetBondTypeAsDouble() >= 2:  # using GetBondTypeAsDouble() is a trick in RDKit
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    def find_glycerol_backbone(mol):
        # Iterate over all sets of three distinct carbons that are connected in a chain: a - b - c.
        carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
        n = len(carbons)
        for i in range(n):
            a = carbons[i]
            for j in range(n):
                if j == i:
                    continue
                b = carbons[j]
                # Check a and b are bonded:
                if not mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx()):
                    continue
                for k in range(n):
                    if k == i or k == j:
                        continue
                    c = carbons[k]
                    # Check b and c are bonded (chain: a - b - c)
                    if not mol.GetBondBetweenAtoms(b.GetIdx(), c.GetIdx()):
                        continue
                    # Now we have a chain a - b - c.
                    # For each of a, b, c, collect oxygen neighbors that are not in {a,b,c}:
                    def get_extras(carbon):
                        extras = []
                        for nbr in carbon.GetNeighbors():
                            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in {a.GetIdx(), b.GetIdx(), c.GetIdx()}:
                                extras.append(nbr)
                        return extras
                    extras_a = get_extras(a)
                    extras_b = get_extras(b)
                    extras_c = get_extras(c)
                    # We require that across these three carbons we have exactly 3 oxygen substituents,
                    # and among them exactly one is phosphate (attached to P) and exactly two are fatty ester oxygens.
                    oxygens = extras_a + extras_b + extras_c
                    if len(oxygens) != 3:
                        continue
                    phosphate_count = 0
                    fatty_count = 0
                    for o in oxygens:
                        if is_phosphate_oxygen(o):
                            phosphate_count += 1
                        elif is_fatty_ester_oxygen(o):
                            fatty_count += 1
                    if phosphate_count == 1 and fatty_count == 2:
                        return True
        return False

    if not find_glycerol_backbone(mol):
        return False, "No proper glycerol backbone with 1 phosphate and 2 fatty acyl ester substituents found"
    
    return True, ("Molecule has one phosphorus, no nitrogen, exactly two fatty acid ester groups, "
                  "and a glycerol backbone with the expected substitution pattern – consistent with phosphatidic acid.")

# For testing purposes, one might call the function with one of the provided SMILES.
if __name__ == "__main__":
    # Example test: PA(14:0/16:1(9Z))
    test_smiles = "[C@](COC(=O)CCCCCCCCCCCCC)(OC(=O)CCCCCCC/C=C\\CCCCCC)([H])COP(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)