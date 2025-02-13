"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA that results from the formal condensation 
            of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.
            
A molecule in this class must contain:
 - A recognizable CoA moiety (here, a fragment typical in many CoA derivatives)
 - A thioester group (C(=O)S) linking an acyl chain to CoA.
 - In the acyl chain attached via the thioester, when numbering the carbonyl carbon as C1,
   the third carbon (C3) – that is, the β-carbon (two bonds away from the carbonyl carbon) – 
   must carry a hydroxyl (-OH) substituent.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES string.
      2. Checks for a Coenzyme A fragment (a typical substructure present in CoA derivatives)
         by matching a characteristic SMARTS.
      3. Searches for thioester groups (C(=O)S) in the molecule.
      4. For each thioester found, traverses the acyl chain: finds the alpha-carbon (the one 
         attached to the carbonyl carbon excluding the carbonyl oxygen and the sulfur) and then 
         its carbon neighbor (the beta carbon). The beta carbon in a 3-hydroxy fatty acyl should have an -OH.
      5. If such a pattern is found, the molecule is classified as part of the 3-hydroxy fatty acyl-CoA class.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as 3-hydroxy fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    
    # Try to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for Coenzyme A substructure.
    # Many CoA derivatives include a fragment like "SCCNC(=O)CCNC(=O)" part of the pantetheine unit.
    coa_smarts = "[S]CCNC(=O)CCNC(=O)"
    coa_frag = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_frag):
        return False, "Coenzyme A fragment missing"
    
    # Look for thioester groups defined as "C(=O)S".
    thioester_smarts = "C(=O)S"
    thioester_frag = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_frag)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found"
    
    # For each thioester, check if the acyl chain has a hydroxyl on the beta carbon (i.e. C3 of fatty acyl)
    found_3hydroxy = False
    for match in thioester_matches:
        # match returns a tuple of indices for atoms matching the pattern:
        # match[0]: the carbonyl carbon (C1), match[1]: the carbonyl oxygen, match[2]: the sulfur.
        c_co = mol.GetAtomWithIdx(match[0])
        
        # Identify the neighbor of the carbonyl carbon that is on the acyl chain.
        # Exclude the oxygen (from the C=O double bond) and the sulfur.
        acyl_neighbors = [nbr for nbr in c_co.GetNeighbors() if nbr.GetIdx() not in match[1:]]
        if not acyl_neighbors:
            continue
        # Choose one neighbor that is a carbon (this should be the α‐carbon, C2)
        alpha = None
        for nbr in acyl_neighbors:
            if nbr.GetAtomicNum() == 6:
                alpha = nbr
                break
        if alpha is None:
            continue
        
        # From the alpha-carbon, get the next carbon along the acyl chain (β‐carbon, C3) 
        alpha_neighbors = [nbr for nbr in alpha.GetNeighbors() if nbr.GetIdx() != c_co.GetIdx()]
        if not alpha_neighbors:
            continue
        beta = None
        for nbr in alpha_neighbors:
            if nbr.GetAtomicNum() == 6:
                beta = nbr
                break
        if beta is None:
            continue
        
        # Check if the beta carbon has a hydroxyl (-OH) substituent.
        hydroxyl_found = False
        for nbr in beta.GetNeighbors():
            # We expect an -OH group: oxygen, typically with only one connection (the hydrogen is implicit)
            if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                hydroxyl_found = True
                break
        
        if hydroxyl_found:
            found_3hydroxy = True
            break  # We only need one valid acyl chain
        
    if not found_3hydroxy:
        return False, "Thioester found but no acyl chain with a hydroxyl on the beta carbon (C3)"
    
    return True, "Molecule contains a CoA moiety and a thioester group with a 3-hydroxy acyl chain"