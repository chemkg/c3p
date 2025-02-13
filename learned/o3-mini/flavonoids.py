"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: Flavonoids (a superclass comprising entities related to flavonoid, isoflavonoid, chalcones, etc.)
Definition (approximate): Organic molecules based on derivatives of a phenyl‐substituted 1-phenylpropane possessing a C15 or C16 skeleton
or structures condensed with C6‐C3 lignan precursors.
Note: This classification is heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid (or closely related) based on its SMILES string.
    
    This function uses two complementary criteria:
      1. A search for a typical flavonoid substructure. We look for either a benzopyran (chroman) motif 
         or for a chalcone-like fragment. These patterns approximate many flavonoid cores.
      2. An analysis of the molecular scaffold (using Murcko scaffolds) to see if it roughly comprises
         15 or 16 carbon atoms (i.e. the C6–C3–C6 core), which is typical for the flavonoid skeleton.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches criteria for flavonoids, False otherwise.
        str: An explanation of the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least two aromatic rings (flavonoids normally have several fused/aromatic rings)
    ri = mol.GetRingInfo()
    aromatic_rings = 0
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings += 1
    if aromatic_rings < 2:
        return False, f"Found only {aromatic_rings} aromatic rings; flavonoids usually contain multiple aromatic rings"
    
    # Define substructure patterns for a benzopyran (chroman) and for a chalcone fragment.
    # The benzopyran pattern here looks for a benzene ring fused to an oxygen-containing heterocycle.
    benzopyran_smarts = "c1ccc2occc2c1"
    benzopyran_pattern = Chem.MolFromSmarts(benzopyran_smarts)
    
    # The chalcone-like pattern is an approximate pattern for two aromatic rings connected by a carbonyl and a short chain.
    chalcone_smarts = "c1ccc(cc1)C(=O)CC=c2ccccc2"
    chalcone_pattern = Chem.MolFromSmarts(chalcone_smarts)
    
    # Check if either pattern is found in the molecule.
    has_flavonoid_motif = mol.HasSubstructMatch(benzopyran_pattern) or mol.HasSubstructMatch(chalcone_pattern)
    
    # Get the Murcko scaffold to remove extra substituents.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error computing Murcko scaffold: {str(e)}"
    
    # Count the number of carbon atoms in the scaffold.
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if the scaffold is in the right ballpark (C15 or C16 core)
    has_correct_scaffold = scaffold_carbons in (15, 16)
    
    # Provide reasoning based on the two analyses.
    if has_flavonoid_motif:
        reason = "Molecule contains a flavonoid motif (benzopyran or chalcone-like fragment)."
    elif has_correct_scaffold:
        reason = f"Molecule has a Murcko scaffold with {scaffold_carbons} carbon atoms, consistent with a C15/C16 flavonoid core."
    else:
        # If neither test gives a positive clue then reject.
        reason = ("No clear flavonoid features found: neither a benzopyran/chalcone motif was detected "
                  f"nor was a core with 15 or 16 carbons identified (scaffold carbons = {scaffold_carbons}).")
        return False, reason

    # Optionally, additional checks – for example, the molecule should be organic.
    # We can check that most heavy atoms are C, O, N or common halogens.
    allowed_atomic_nums = {6, 7, 8, 9, 17, 35}  # C, N, O, F, Cl, Br
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 20 and atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()} that is unusual in flavonoids"

    return True, f"Flavonoid-like: {reason}"

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES for a flavonoid (azaleatin, for instance)
    test_smiles = "COc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O)c(=O)c12"
    result, message = is_flavonoids(test_smiles)
    print(result, message)