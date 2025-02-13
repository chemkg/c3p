"""
Classifies: CHEBI:36916 cation
"""
#!/usr/bin/env python
"""
Classifies: cation
Definition:
  A monoatomic or polyatomic species having one or more elementary charges of the proton.
  
Improved criteria:
  1. The molecule must parse to a single fragment.
  2. Very small molecules (fewer than 4 heavy atoms) are rejected unless they are known metal ions.
  3. If the overall formal charge is positive, the molecule is accepted.
  4. If the overall charge is zero, then (a) if a positive center exists that is not directly bonded to any negative
     AND the balance of charged atoms is not exactly 1:1, we regard it as a cation, OR
     (b) if the molecule exactly balances positive and negative charges, we check for a known phosphocholine moiety.
     
Note: This is a heuristic classifier. Some borderline cases (including zwitterions) may be mis‐classified.
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    
    The function works as follows:
      1. Parse the SMILES string into an RDKit molecule.
      2. If the molecule consists of more than one disconnected fragment, it is not considered a single ionic species.
      3. If the molecule is very small (fewer than 4 heavy atoms) and is not a known metal ion, it is rejected.
      4. The overall net formal charge is computed. A net positive charge qualifies immediately.
      5. If the net charge is zero then two routes are taken:
          (a) If positive atoms outnumber negative ones and at least one positively charged atom has no 
              immediate negative neighbor, then the molecule is accepted.
          (b) Otherwise, if the positive and negative counts are equal, we look for a characteristic 
              phosphocholine substructure (common in some lipids known to behave as isolated cations).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple of True/False and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Check if molecule is a single fragment.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) != 1:
        return False, "Multiple disconnected fragments found (likely a salt pair); cannot classify as a cation."
    
    # Reject very small molecules unless they are known metal ions.
    # (For example, "[Co+]" should be accepted but O=[I+]=O should not.)
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(heavy_atoms) < 4:
        # Define allowed metal ions (if all heavy atoms are from this set then we allow the species)
        allowed_metals = {"Co", "Cr", "Na", "K", "Li", "Ag", "Au", "Hg"}
        if not all(atom.GetSymbol() in allowed_metals for atom in heavy_atoms):
            return False, f"Too few heavy atoms ({len(heavy_atoms)}) in non-metal species; not classified as cation."
    
    # Calculate overall net formal charge as well as counts of positive and negative centers.
    net_charge = 0
    pos_count = 0
    neg_count = 0
    for atom in mol.GetAtoms():
        fc = atom.GetFormalCharge()
        net_charge += fc
        if fc > 0:
            pos_count += 1
        elif fc < 0:
            neg_count += 1

    # If overall net charge is positive, classify as cation.
    if net_charge > 0:
        return True, f"Overall net formal charge is {net_charge}."
    
    # If overall net charge is negative, it is not a cation.
    if net_charge < 0:
        return False, f"Overall net formal charge is {net_charge} (anionic species)."
    
    # Here net_charge == 0.
    # Option 1: If there are more positive than negative atoms, check if any positive atom is "unpaired" (neighbors have no negative).
    if pos_count > neg_count:
        for atom in mol.GetAtoms():
            if atom.GetFormalCharge() > 0:
                # Check if any neighbor is negative.
                if all(nb.GetFormalCharge() >= 0 for nb in atom.GetNeighbors()):
                    return True, (f"Neutral overall but has unpaired positive center: {atom.GetSymbol()}({atom.GetFormalCharge()}). "
                                  f"Positives: {pos_count}, Negatives: {neg_count}.")
        # If none of the positives are isolated, then do not classify as cation.
        return False, f"Overall net formal charge is 0 with all positive centers paired with negative neighbors."
    
    # Option 2: If positive and negative counts are exactly equal (a balanced zwitterion),
    # look for a specific phosphocholine-type moiety that is known to behave as a cation.
    if pos_count == neg_count and pos_count > 0:
        # Define a SMARTS pattern for a phosphocholine moiety:
        # This pattern looks for a phosphate group bonded to an [O-] and then an ethoxy chain ending in a tetraalkylammonium.
        phosphocholine_smarts = "[O-]P(=O)([O-])OCC[N+](C)(C)C"
        pchol_pattern = Chem.MolFromSmarts(phosphocholine_smarts)
        if mol.HasSubstructMatch(pchol_pattern):
            return True, ("Neutral overall but contains a phosphocholine moiety "
                          "([O-]P(=O)([O-])OCC[N+](C)(C)C) that behaves as an unpaired cationic center.")
    
    return False, f"No unpaired positive centers found. Overall net formal charge is {net_charge} with positives: {pos_count} and negatives: {neg_count}."

# Example usage: (this block can be commented out when importing as a module)
if __name__ == "__main__":
    test_smiles_list = [
        "[H][C@@]12[NH+]3CC[C@]11C(Nc4cc(O)ccc14)=C(C[C@]2(CC)C=CC3)C(=O)OC",  # 16-hydroxytabersoninium (net +)
        "C([C@](CO/C=C\\CCCCCCCCCCCCCCCC)([H])OC(CCC/C=C\\C/C=C/CCCCC)OP([O-])(=O)OCC[N+](C)(C)C",  # PC(P-18:0/20:5(5Z,8Z,11Z,14Z,17Z))
        "CCCCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC)COP([O-])(=O)OCC[N+](C)(C)C",  # 1-[(9Z,12Z)-octadecadienoyl]-2-hexadecanoyl-sn-glycero-3-phosphocholine
        "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCCCCC)([O-])=O",  # PC(19:0/15:1(9Z))
        "CN(C)c1ccc(cc1)C(=C1C=CC(C=C1)=[N+](C)C)c1ccc(cc1)[N+](C)(C)C",  # methyl green(2+)
        "CCN1c2cc3Oc4cc5=[N+](CC)C(C)(C)C=Cc5cc4=C(c3cc2C(C)=CC1(C)C)c1cc(ccc1C(O)=O)C(O)=O",  # ATTO 590 para-isomer(1+)
        "C(CCC)CC[NH3+]",  # hexan-1-aminium  (net +)
        "C1(=C(C=C(C=C1OC)CC[NH3+])O)O",  # 3,4-dihydroxy-5-methoxyphenethylaminium (net +)
        "[NH3+][C@@H]([C@H]1CCNC(=[NH2+])N1)C([O-])=O",  # (2S,3R)-capreomycidine(1+)
        "[Co+]",  # cobalt(1+)
        "O=[I+]=O",  # dioxidoiodine(1+): should be rejected by small‐molecule rule.
        "B#[O+]",  # oxidoboron(1+): should be rejected by small‐molecule rule.
        "[Na+].CCCC([O-])=O",  # salt pair: should be rejected.
    ]
    
    for s in test_smiles_list:
        result, reason = is_cation(s)
        print(f"SMILES: {s}\n  -> Cation: {result}\n  Reason: {reason}\n")