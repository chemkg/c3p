"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: phenylpropanoid compounds
Definition: Any organic aromatic compound with a structure based on a 
phenylpropane (C6–C3) skeleton. This class includes naturally occurring 
phenylpropanoid esters, flavonoids, anthocyanins, coumarins and many small 
phenolic molecules as well as their semi‐synthetic and synthetic analogues.
Phenylpropanoids are also precursors of lignin.

This implementation uses several SMARTS patterns to capture characteristic substructures:
  • Coumarin (and 4‐hydroxycoumarin) 
  • Cinnamic acid derivatives – now requiring a CH₂ (benzylic) group attached to benzene
  • A benzene with an unsaturated “cinnamyl” chain (again, demanding a CH₂ linker)
  • Flavanone, isoflavone, flavonol and flavan‐3‐ol cores

If none of these patterns directly match on the molecule, then the Murcko scaffold (the core ring system) is computed 
and the patterns are re‐applied. This helps rescue many flavonoid glycosides or heavily decorated structures.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.

    The heuristic requires that the molecule has at least one benzene ring and then 
    looks for one of a list of SMARTS motifs characteristic of phenylpropanoid-type 
    molecules:
      - Coumarin structure.
      - Cinnamic acid derivative (requires a benzene ring linked via a CH2 to a vinyl and then a carbonyl).
      - Benzene with an unsaturated (cinnamyl) chain (again via a CH2 linker).
      - Flavanone core.
      - Isoflavone core.
      - Flavonol core.
      - Flavan-3-ol core.
      
    If no direct match is found on the full molecule then the Murcko scaffold is computed and the 
    same patterns are re-applied. This helps catch cases where sugars or extra substituents “hide” the main core.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if classified as a phenylpropanoid, False otherwise.
      str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule to compute aromaticity and valence.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Basic check: at least one benzene ring must be present.
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene ring detected, so not a phenylpropanoid"
    
    # A list of (description, SMARTS) pairs for phenylpropanoid cores.
    patterns = [
        # Coumarin: typical coumarin or 4-hydroxycoumarin core.
        ("Coumarin structure", "O=c1ccc2oc(=O)cc2c1"),
        # Cinnamic acid derivative: a benzene ring connected via a CH2 (benzylic) group to a vinyl and then a carbonyl (acid or anion)
        ("Cinnamic acid derivative", "c1ccccc1[CH2]/C=C/C(=O)[O*]"),
        # Benzene with unsaturated (cinnamyl) chain: requires an explicit CH2 linker.
        ("Benzene with unsaturated (cinnamyl) chain", "c1ccccc1[CH2]/C=C/\\C"),
        # Flavanone core (typical bicyclic flavonoid): note that stereochemistry is not forced.
        ("Flavanone core", "c1cc(c(cc1)O)C2=CC(=O)OC3=CC=CC=C23"),
        # Isoflavone core.
        ("Isoflavone core", "c1ccc2C(=O)oc3ccccc3c2c1"),
        # Flavonol core: a typical 3-hydroxyflavone scaffold.
        ("Flavonol core", "c1cc(c(c(c1O)O)C2=CC(=O)C(O)=C2)O"),
        # Flavan-3-ol core (e.g. catechins).
        ("Flavan-3-ol core", "c1cc(c(c1O)O)C2CC(O)C(O)C2")
    ]
    
    # Try matching each pattern.
    for desc, smarts in patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue  # skip if SMARTS is invalid (should not happen)
        if mol.HasSubstructMatch(patt):
            return True, f"Matches pattern: {desc}"
    
    # If no pattern was found on the full molecule, try using the Murcko scaffold – 
    # this can remove peripheral sugars or other substituents.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        # Sanitize the scaffold if not already done.
        Chem.SanitizeMol(scaffold)
    except Exception:
        scaffold = None
    if scaffold is not None:
        for desc, smarts in patterns:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            if scaffold.HasSubstructMatch(patt):
                return True, f"Matches pattern on scaffold: {desc}"
    
    # If nothing was found, the molecule is not classified as a phenylpropanoid.
    return False, "No phenylpropanoid-characteristic substructure found"

# Example usage (for testing—remove or comment out if used as a module):
if __name__ == "__main__":
    test_smiles = [
        # Compounds expected to be phenylpropanoids (from the provided examples):
        "C1(C2=C(O[C@@H](C1)C3=C(C=C(C=C3)O)O)C(=C(C=C2O)O)C[C@@H](CCC(=C)C)C(=C)C)=O",  # remangiflavanone B
        "COc1ccc(cc1OC)[C@H]1O[C@H]([C@H](C)[C@@H]1C)c1ccc2OCOc2c1",                   # futokadsurin C
        "COc1cc2C[C@H]3COC(=O)[C@@H]3Cc3c(OC)c(OC)c(OC)cc3-c2cc1OC",                   # Neoisostegane
        "COC1=C(O)C=CC(=C1)[C@H]1OC2=CC(O)=CC(O)=C2C(=O)[C@@H]1O",                       # (+)-dihydroisorhamnetin
        "O1C(C(O)C(O)C(OC(=O)/C=C/C2=CC=C(O)C=C2)C1OC=3C=C4OC=C(C(=O)C4=C(O)C3)C5=CC=C(O)C=C5)CO",  # Genistein 7-O-(2-p-coumaroylglucoside)
        # Compounds that should NOT be classified as phenylpropanoids:
        "OC(=Cc1ccccc1)C([O-])=O",  # enol-phenylpyruvate (previous false positive)
        "CCO"  # ethanol
    ]
    
    for sm in test_smiles:
        result, reason = is_phenylpropanoid(sm)
        print(f"SMILES: {sm}\nClassified as phenylpropanoid? {result}\nReason: {reason}\n")