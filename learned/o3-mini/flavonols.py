"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: Flavonols (CHEBI: flavonols)
Definition: Any hydroxyflavone in which the ring hydrogen at position 3 of the heterocyclic (chromen) ring is replaced by a hydroxy group.
Our strategy:
  - Construct a SMARTS pattern corresponding to the core structure of 3-hydroxyflavone,
    i.e. a chromen-4-one skeleton with a hydroxy group at the position that would be 3 in flavonols.
  - Do not force the B-ring to be an unsubstituted phenyl; allow arbitrary substitution.
  - Remove stereochemistry from the query molecule to improve matching.
"""

from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol (i.e. contains a 3-hydroxyflavone scaffold)
    based on its SMILES string. A flavonol is defined as any hydroxyflavone in which the ring 
    hydrogen at position 3 is replaced by a hydroxy group.
    
    Our strategy is to search for the 3-hydroxyflavone core using a SMARTS pattern.
    To improve matching for derivatives, stereochemical labels are removed and useChirality is set False.
    
    The SMARTS we use is based on the chromen-4-one (flavone) nucleus with an -OH at position 3.
    We relax the substitution pattern on the B-ring.
    
    For this pattern:
      - "O=c1c(O)cc2oc(c1)cc2" 
        It finds:
          * a carbonyl group (O=) on an aromatic ring (position 4)
          * a fused heterocycle (the chromen system)
          * a hydroxyl group (O) directly attached to the ring carbon that, in a typical flavonol, is at position 3.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavonol, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemical information to avoid mismatches due to chiral centers.
    Chem.RemoveStereochemistry(mol)
    
    # Define a SMARTS pattern for the 3-hydroxyflavone core.
    # We drop the explicit specification of the phenyl (B-ring) to allow for derivatives.
    flavonol_core_smarts = "O=c1c(O)cc2oc(c1)cc2"
    flavonol_pattern = Chem.MolFromSmarts(flavonol_core_smarts)
    if flavonol_pattern is None:
        return None, None  # In case the SMARTS construction fails.
    
    # Perform the substructure search ignoring chirality.
    if mol.GetSubstructMatches(flavonol_pattern, useChirality=False):
        return True, "Molecule contains the 3-hydroxyflavone core."
    else:
        return False, "Molecule does not contain the required 3-hydroxyflavone (flavonol) core."

# For testing purposes:
if __name__ == "__main__":
    # List of test SMILES strings from the provided examples.
    test_smiles = [
        "OC[C@H]1O[C@@H](Oc2ccc(cc2)-c2oc3cc(O)cc(O)c3c(=O)c2O",  # kaempferol 4'-O-beta-D-glucopyranoside
        "C12=C(OC(C3=CC=C(OC)C=C3)=C(C1=O)O)C=C(OC)C=C2O",            # 7,4'-dimethylkaempferol
        "COc1cc(O)c2c(c1)oc(-c1ccc(OC)c(O)c1)c(O)c2=O",              # ombuin
        "c1ccccc1-c2oc3ccccc3c(=O)c2",                                # A flavone lacking 3-OH (should be False)
        "COc1cc(O)c2c(c1)oc(-c1ccccc1)c(O)c2=O"                       # kaempferide (as example)
    ]
    
    for smi in test_smiles:
        result, reason = is_flavonols(smi)
        print(f"SMILES: {smi}")
        print(f"Flavonol? {result} ({reason})\n")