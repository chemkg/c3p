"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione)
skeleton with a substituent attached on the 10 position (i.e. a flavin).
This approach uses a two‐step strategy:
  1. A lenient SMARTS is used to find an isoalloxazine core.
  2. The candidate N atom (mapped as N10 in the query) is inspected:
       • In an unsubstituted isoalloxazine this N is drawn as [nH] (with one hydrogen).
         In a flavin the substituent replaces that hydrogen so that the N has 0 hydrogens.
       • We further require that a substituent attached at that N outside of the core is carbon‐based.
Note: This is a heuristic and may not guarantee 100% accuracy.
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin derivative (i.e. a derivative of the dimethylisoalloxazine
    skeleton with a substituent at the 10 position) based on its SMILES string.
    
    Our strategy is:
      1. Parse the SMILES.
      2. Look for the isoalloxazine core using a lenient SMARTS.
         In an unsubstituted isoalloxazine the N at the 10-position is drawn as [nH].
         But as a flavin derivative that N loses its hydrogen when substituted.
         To allow matching both cases we write our SMARTS as:
           "n1c2cc(cc2)nc2c1nc(=O)[n:10]c2=O"
         (i.e. we merely tag the N expected at the 10 position without imposing a hydrogen).
      3. For the matched core, examine the candidate N (mapped as 10):
           - Check the total number of attached hydrogens (both explicit and implicit).
             In a flavin derivative the N should have 0 hydrogens (the substituent has replaced it).
           - Verify that at least one neighbor of that N, which is not part of the core match,
             is a carbon (atomic number 6) as expected if a ribityl or similar chain is present.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a flavin derivative, False otherwise.
       str: A reason describing the decision.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a lenient SMARTS pattern to capture the isoalloxazine core.
    # The pattern is modified relative to the previous version to not require [nH]
    # at the 10-position. We tag that nitrogen as [n:10].
    core_smarts = "n1c2cc(cc2)nc2c1nc(=O)[n:10]c2=O"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in the SMARTS pattern for the flavin core."
    
    # Find substructure matches of the isoalloxazine core in the molecule.
    matches = mol.GetSubstructMatches(core_query)
    if not matches:
        return False, "Isoalloxazine core not found."
    
    # For each match, check the condition on N10.
    for match in matches:
        # Find the atom in the query with mapping number 10.
        n10_query_idx = None
        for atom in core_query.GetAtoms():
            if atom.HasProp("molAtomMapNumber") and atom.GetProp("molAtomMapNumber") == "10":
                n10_query_idx = atom.GetIdx()
                break
        if n10_query_idx is None:
            return False, "The query pattern does not contain the expected atom mapping for N10."
        
        # Get the corresponding atom in the molecule.
        n10_mol_idx = match[n10_query_idx]
        n10_atom = mol.GetAtomWithIdx(n10_mol_idx)
    
        # In an unsubstituted isoalloxazine, the N10 would be drawn with an attached hydrogen ([nH]).
        # In a flavin derivative the substituent at N10 replaces that hydrogen, so total H count should be 0.
        if n10_atom.GetTotalNumHs() > 0:
            return False, ("Isoalloxazine core found but the N10 atom still has a hydrogen; "
                           "expected substituent replacing the hydrogen for a flavin.")
    
        # Retrieve the set of core atom indices from our current match.
        core_atom_indices = set(match)
        # Now verify that at least one neighbor (outside of the core) is carbon-based.
        valid_substituent = False
        for neighbor in n10_atom.GetNeighbors():
            if neighbor.GetIdx() not in core_atom_indices:
                if neighbor.GetAtomicNum() == 6:
                    valid_substituent = True
                    break
        if not valid_substituent:
            return False, "Substituent at N10 is not carbon‐based as expected for a flavin."
    
        return True, "Molecule contains the expected isoalloxazine core with a substituent on N10."
    
    return False, "No valid flavin substructure match found."

# Example usage:
if __name__ == "__main__":
    # Testing with several example SMILES strings.
    test_smiles = {
        "8-formyl-8-demethylriboflavin 5'-phosphate":
            "C(N1C=2C(=NC3=C1C=C(C(=C3)C)C(=O)[H])C(NC(N2)=O)=O)[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)O",
        "Riboflavol":
            "OC(CN1C2=C(N=C3C1=NC(=O)N(O)C3=O)C=C(C(=C2)C)C)C(O)C(O)CO",
        "riboflavin cyclic 4',5'-phosphate":
            "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H](O)[C@H](O)[C@H]3COP(O)(=O)O3)c2cc1C",
        "prenyl-FMN":
            "C=12N=C(NC(C1[N+]=3C=4C(N2C[C@@H]([C@@H]([C@@H](COP(O)(O)=O)O)O)O)=CC(=C(C4C(C)(CC3)C)C)C)=O)[O-]",
        "FADH(.)":
            "Cc1cc2[N]c3c([nH]c(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)OP(O)(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C"
    }
    for name, smi in test_smiles.items():
        flag, reason = is_flavin(smi)
        print(f"{name}: {flag} – {reason}")