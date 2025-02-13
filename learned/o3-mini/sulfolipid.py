"""
Classifies: CHEBI:61384 sulfolipid
"""
#!/usr/bin/env python
"""
Classifies: sulfolipid
Definition: A compound containing a sulfonic acid residue joined (via an appropriate linkage)
            to a lipid (i.e., having a long aliphatic chain).
            
This function uses RDKit to verify that the molecule contains a sulfo group (–OS(=O)(=O)O or
its deprotonated variant) and a long chain (a substructure match of 10 consecutive carbon atoms).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if the given SMILES string corresponds to a sulfolipid.
    
    The strategy is:
      1. Parse the SMILES.
      2. Look for a sulfonic acid (sulfo) moiety. Here we look for either the protonated
         form "OS(=O)(=O)O" or the deprotonated form "OS(=O)(=O)[O-]".
      3. Look for a lipid chain. We require a long contiguous carbon chain (here 10 carbons in sequence).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a sulfolipid, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Define SMARTS pattern(s) for sulfo group.
    # Many sulfolipids show a sulfo group as –OS(=O)(=O)O or –OS(=O)(=O)[O-]
    sulfo_protonated = Chem.MolFromSmarts("OS(=O)(=O)O")
    sulfo_deprotonated = Chem.MolFromSmarts("OS(=O)(=O)[O-]")
    
    has_sulfo = False
    if mol.HasSubstructMatch(sulfo_protonated):
        has_sulfo = True
    elif mol.HasSubstructMatch(sulfo_deprotonated):
        has_sulfo = True
        
    if not has_sulfo:
        return False, "No sulfonic acid (sulfo) group found."
    
    # Define a pattern for a long aliphatic chain.
    # Here we use a simple SMARTS for 10 consecutive aliphatic carbons.
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long aliphatic (lipid) chain found."
    
    # Optionally, you might want to check some overall molecular properties (such as molecular weight,
    # rotatable bonds etc.) to further support the lipid character. For now, we only apply the two substructure
    # searches.
    
    return True, "Contains a sulfonic acid (sulfo) moiety and a long aliphatic (lipid) chain."

# The module can be tested by calling is_sulfolipid with known sulfolipid SMILES.
# For example:
# result, reason = is_sulfolipid("CCCCCCCCCCCCCCCC(=O)O[C@H]1...")  # (example SMILES)
# print(result, reason)