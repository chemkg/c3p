"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine
Definition:
  A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
  That is, a glycerol backbone in which the sn-3 hydroxyl is phosphorylated and further esterified with L-serine,
  and the sn-1 and sn-2 hydroxyls are esterified with fatty acyl chains.
  
This program uses RDKit and two highly specific SMARTS patterns (allowing for optionally deprotonated groups)
to recognize the head‐group of a 3‑sn‐phosphatidyl‐L‐serine.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    This function verifies that the molecule contains a glycerophosphoserine head group in which:
      - the glycerol backbone has three –OH groups,
      - the sn-3 hydroxyl is replaced by a phosphate that is further O‑linked to an L‑serine unit,
      - and the sn-1 and sn-2 hydroxyls are esterified with fatty acyl chains.
    
    In this implementation we use two complementary (mirror image) SMARTS patterns.
    These patterns are written to (a) allow common representations such as deprotonated oxygens
    and (b) “lock in” the connectivity of the glycerol–phosphate–serine fragment with two ester linkages.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as a 3-sn-phosphatidyl-L-serine, False otherwise.
      str: Reason for the classification.
    """
    # Attempt to parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (Optional) use a molecular weight filter to rule out very small molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low to be a phosphatidylserine"
    
    # Here we define two SMARTS patterns for the entire head-group.
    # The intended fragment is that of a glycerol backbone in which:
    #   * The sn-1 –OH is esterified: shown as O C(=O)[*]  (acyl chain; [*] can be any chain)
    #   * The sn-3 –OH is replaced by a phosphate which (allowing optionally for deprotonation)
    #     is further O‑linked to an L‐serine moiety (with a chiral center and an amino and carboxyl group).
    # We allow alternatives for protonated/deprotonated oxygens by using [$([O]),$([O-])].
    #
    # Pattern explanation (ps_full_smarts1):
    #   O C[C@H](    -- the glycerol backbone central carbon with defined chirality
    #        COP(=O)([$([O]),$([O-])])([$([O]),$([O-])])
    #               C[C@H](N)C(=O)[$([O]),$([O-])]
    #         )
    #         OC(=O)[*]  -- the ester at the remaining glycerol –OH (acyl chain)
    #
    # The second pattern is the mirror image with [C@@H] at the glycerol center.
    ps_full_smarts1 = ("OC[C@H](COP(=O)([$([O]),$([O-])])([$([O]),$([O-])])"
                        "C[C@H](N)C(=O)[$([O]),$([O-])])OC(=O)[*]")
    ps_full_smarts2 = ("OC[C@@H](COP(=O)([$([O]),$([O-])])([$([O]),$([O-])])"
                        "C[C@H](N)C(=O)[$([O]),$([O-])])OC(=O)[*]")
    
    # Build the query molecules from the SMARTS
    query1 = Chem.MolFromSmarts(ps_full_smarts1)
    query2 = Chem.MolFromSmarts(ps_full_smarts2)
    
    if query1 is None and query2 is None:
        # In the unlikely event that our SMARTS failed to compile,
        return False, "Error in head-group SMARTS definitions"
    
    # Check if at least one of the patterns is found in the molecule
    if mol.HasSubstructMatch(query1) or mol.HasSubstructMatch(query2):
        return True, ("Molecule contains a glycerophosphoserine head group with two acyl ester substituents "
                      "at the glycerol sn-1 and sn-2 positions")
    else:
        return False, "Phosphatidyl-L-serine head group pattern not found"

# Optional testing:
# test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCC"
# print(is_3_sn_phosphatidyl_L_serine(test_smiles))