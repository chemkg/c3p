"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid
Definition: Any member of class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage to a carbohydrate part (usually a mono-, di- or tri-saccharide). Some substances classified as bacterial glycolipids have the sugar part acylated by one or more fatty acids and the glycerol part may be absent.
Heuristics used:
 - Look for a carbohydrate moiety via a pyranose-like ring pattern.
 - Look for at least two ester groups (as a proxy for acyl chains on a glycerol backbone)
 - Alternatively, allow the possibility of a glycosphingolipid or bacterial glycolipid by checking for a long aliphatic chain (heuristically, the substring "CCCCCCCC" in the SMILES) or an amide bond with long chain.
Note: This is an approximate classification.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    
    A glycolipid (by our working definition) is a compound that contains:
     - a carbohydrate (sugar) moiety joined glycosidically (via an O-linkage) and 
     - a lipid part, usually in the form of a 1,2-di-O-acylglycerol (two ester groups),
       or in other cases a sphingolipid type (amide + long alkyl chain) or a long aliphatic chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our glycolipid criteria, False otherwise.
        str: Reason for the classification decision.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Heuristic 1: Detect a carbohydrate moiety ---
    # Here we use a SMARTS pattern for a pyranose-like ring.
    # This pattern looks for a 6-membered ring containing an oxygen with several OH groups.
    sugar_pattern = Chem.MolFromSmarts("OC1OC(O)C(O)C(O)C1O")
    sugar_found = mol.HasSubstructMatch(sugar_pattern)
    
    if not sugar_found:
        return False, "No clear carbohydrate (sugar) moiety detected"
    
    # --- Heuristic 2: Look for acyl (ester) groups ---
    # A simple ester group pattern is defined as a carbonyl directly connected to an oxygen.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    n_esters = len(ester_matches)
    
    # Also look for amide bonds as one might find in glycosphingolipids.
    amide_pattern = Chem.MolFromSmarts("N[C;!$(C=O)]C(=O)")
    amide_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("NC(=O)"))
    n_amides = len(amide_matches)
    
    # --- Heuristic 3: Look for a long aliphatic chain ---
    # We use a simple substring search as a crude indicator.
    long_chain = "CCCCCCCC" in smiles  # at least 8 contiguous carbons
    
    # Decision logic:
    # Typical glycolipids (1,2-di-O-acylglycerols) should have at least 2 ester groups.
    # However, some glycolipids (e.g. glycosphingolipids or bacterial glycolipids) lack a typical glycerol and may show an amide bond or a long chain.
    if n_esters >= 2:
        return True, f"Contains carbohydrate moiety and {n_esters} ester group(s) (di-acyl structure expected)"
    elif n_amides >= 1 and long_chain:
        return True, f"Contains carbohydrate moiety and an amide bond with a long acyl chain"
    elif long_chain:
        # Even if ester/amide count is low, a long chain and sugar might indicate a glycolipid.
        return True, "Contains carbohydrate moiety and a long aliphatic chain, possible glycolipid variant"
    else:
        return False, "Carbohydrate is present but no clear lipid (diacyl or long chain) feature detected"

# Example usage:
if __name__ == '__main__':
    # Some test SMILES strings from the examples (only one or two provided here for demonstration)
    test_smiles = "OC1OC(O)C(O)C(O)C1OCC(=O)CCCCCCCCCCCC"  # This is a very simplified glycolipid-like structure.
    classification, reason = is_glycolipid(test_smiles)
    print("Is glycolipid?", classification)
    print("Reason:", reason)