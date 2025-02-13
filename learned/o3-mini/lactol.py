"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol
Definition:
  Cyclic hemiacetals formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.
  They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.
  
This program uses a heuristic substructure search. In addition to a SMARTS pattern that defines a lactol-like 
motif (a ring sp3-hybridized carbon not a carbonyl center, bonded to an exocyclic hydroxyl group and a ring oxygen),
we add extra filtering:
  - We add explicit hydrogens so that the [OX2H] pattern is reliably matched.
  - For each match, the exocyclic hydroxyl oxygen is examined; if it is “branched” (i.e. attached to another heavy atom 
    besides the lactol carbon) then it is assumed to be part of an O-glycosidic linkage rather than a free hemiacetal.
  - We also require that the lactol carbon and the adjacent ring oxygen belong to a ring of size 5 or 6.
  
Note: This heuristic will not be perfect.
"""

from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol (cyclic hemiacetal) based on its SMILES string.
    A lactol results from the intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl
    and appears as a 1-oxacycloalkan-2-ol (or its unsaturated analogue). 
    
    The approach uses:
      - Adding explicit hydrogens (to capture the hydroxyl hydrogen).
      - A SMARTS pattern that looks for a ring carbon (sp3 and not a carbonyl) that is attached to:
            (a) an exocyclic hydroxyl ([OX2H] not in a ring) and 
            (b) an oxygen that is part of a ring.
      - Post-filtering: the exocyclic hydroxyl should be “free” (i.e. not linked to additional heavy atoms),
            and the lactol carbon and ring oxygen should be together in a 5- or 6-membered ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lactol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit H atoms so that the hydroxyl groups can be recognized
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for a lactol-like center:
    #   - [C;R;X4;!$(C=O)]: A ring carbon (sp3) which is not a carbonyl.
    #   - ([OX2H;!R]): an exocyclic hydroxyl group (oxygen with H, not in a ring)
    #   - [O;R]: a ring oxygen (part of the same cyclic acetal structure).
    lactol_smarts = "[C;R;X4;!$(C=O)]([OX2H;!R])[O;R]"
    lactol_pat = Chem.MolFromSmarts(lactol_smarts)
    if lactol_pat is None:
        return False, "Error in lactol SMARTS pattern"
    
    matches = mol.GetSubstructMatches(lactol_pat)
    if not matches:
        return False, "No lactol moiety pattern (cyclic hemiacetal) detected in the molecule"
    
    # Loop over all matches and apply further checks. Each match is a tuple of atom indices:
    #   match[0]: the candidate lactol carbon,
    #   match[1]: the exocyclic hydroxyl oxygen,
    #   match[2]: the ring oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    for match in matches:
        lactol_c_idx, exo_oh_idx, ring_o_idx = match
        lactol_c = mol.GetAtomWithIdx(lactol_c_idx)
        exo_oh = mol.GetAtomWithIdx(exo_oh_idx)
        ring_o = mol.GetAtomWithIdx(ring_o_idx)
        
        # Check: the exocyclic hydroxyl oxygen should not be further substituted.
        # (In a free hemiacetal, the O should have only the lactol carbon as heavy-atom neighbor.)
        if exo_oh.GetDegree() > 1:
            # Likely part of a glycosidic bond; skip this match.
            continue
        
        # Determine if lactol carbon and ring oxygen belong to a ring of appropriate size (5 or 6 members)
        valid_ring = False
        for ring in ring_info:
            if lactol_c_idx in ring and ring_o_idx in ring:
                if len(ring) in [5, 6]:
                    valid_ring = True
                    break
        
        if valid_ring:
            return True, "Found lactol moiety pattern (cyclic hemiacetal) in the molecule"
    
    return False, "No suitable lactol moiety pattern detected in the molecule"

# Example usage (testing on a couple of SMILES):
if __name__ == "__main__":
    test_smiles = [
        "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha-D-fructopyranose (cyclic hemiacetal)
        "O=C1OC(CO)C(O)C1"  # A small ring model (may or may not be a lactol)
    ]
    for s in test_smiles:
        result, reason = is_lactol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")