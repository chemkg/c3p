"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione)
skeleton with a substituent attached on the 10 position (i.e. a flavin derivative).
This approach uses a two‐step strategy:
  1. A lenient SMARTS is used to find an isoalloxazine core.
  2. For the matched core the candidate N atom (mapped as N10) is verified to be substituted:
       - It must have zero total hydrogens (i.e. the substituent has replaced the expected [nH] proton).
       - At least one neighbor (outside the core match) must be carbon‐based.
Note: This is a heuristic method and may not guarantee 100% accuracy.
"""

from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin derivative based on its SMILES string.
    A flavin derivative is defined as a derivative of the dimethylisoalloxazine 
    (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a carbon substituent 
    attached at the 10 position (replacing the hydrogen normally found on that nitrogen).
    
    Strategy:
      1. Parse the SMILES string.
      2. Look for the isoalloxazine core using a lenient SMARTS pattern.
         The SMARTS pattern is designed to allow the atom at position 10 to be drawn as either 
         aromatic or non-aromatic. (We use "[#7:10]" to indicate any nitrogen atom at that position.)
      3. For any match, check that the candidate N (mapped as atom 10) has no hydrogen attached
         and that at least one of its neighbors outside the matched core is carbon-based.
         
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a flavin derivative, False otherwise.
       str: A reason describing the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a lenient SMARTS pattern for the isoalloxazine core.
    # We allow the N at the 10 position to appear as any nitrogen (aromatic or not)
    # and we tag it with [#7:10].
    # The pattern below loosely captures the fused ring system with two carbonyl groups.
    core_smarts = "n1c2cc(cc2)nc2c1nc(=O)[#7:10]c2=O"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in the SMARTS pattern for the isoalloxazine core."
    
    # Find all substructure matches in the molecule.
    matches = mol.GetSubstructMatches(core_query)
    if not matches:
        return False, "Isoalloxazine core not found."
    
    # Identify the query atom index corresponding to the tagged N at position 10.
    # (We expect the SMARTS to assign mapping number 10 to that atom.)
    n10_query_idx = None
    for atom in core_query.GetAtoms():
        if atom.HasProp("molAtomMapNumber") and atom.GetProp("molAtomMapNumber") == "10":
            n10_query_idx = atom.GetIdx()
            break
    if n10_query_idx is None:
        return False, "SMARTS pattern missing the expected mapping for N10."
    
    # Now evaluate each match.
    for match in matches:
        # Get the corresponding atom from the molecule for the N10 query atom.
        n10_mol_idx = match[n10_query_idx]
        n10_atom = mol.GetAtomWithIdx(n10_mol_idx)
        
        # For a flavin derivative, the N10 should have no hydrogens (i.e. it is substituted).
        if n10_atom.GetTotalNumHs() > 0:
            return False, ("Isoalloxazine core found but N10 still has a hydrogen; "
                           "expected substituent replacing that hydrogen for a flavin derivative.")
        
        # Get the set of atom indices that are part of the matched isoalloxazine core.
        core_atom_indices = set(match)
        
        # Check that at least one neighbor of N10 that is outside the core is carbon based.
        found_carbon_substituent = False
        for neighbor in n10_atom.GetNeighbors():
            if neighbor.GetIdx() not in core_atom_indices and neighbor.GetAtomicNum() == 6:
                found_carbon_substituent = True
                break
        if not found_carbon_substituent:
            return False, "Substituent at N10 is not carbon-based as expected for a flavin derivative."
        
        # If we have come this far, we have a valid match.
        return True, "Molecule contains an isoalloxazine core with a proper substituent at N10."
    
    return False, "No valid flavin substructure match found."

# For testing purposes:
if __name__ == "__main__":
    # Example SMILES strings from the provided list.
    test_smiles = {
        "8-formyl-8-demethylriboflavin 5'-phosphate": "C(N1C=2C(=NC3=C1C=C(C(=C3)C)C(=O)[H])C(NC(N2)=O)=O)[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)O",
        "Riboflavol": "OC(CN1C2=C(N=C3C1=NC(=O)N(O)C3=O)C=C(C(=C2)C)C)C(O)C(O)CO",
        "riboflavin cyclic 4',5'-phosphate": "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H](O)[C@H](O)[C@H]3COP(O)(=O)O3)c2cc1C",
        "prenyl-FMN": "C=12N=C(NC(C1[N+]=3C=4C(N2C[C@@H]([C@@H]([C@@H](COP(O)(O)=O)O)O)O)=CC(=C(C4C(C)(CC3)C)C)C)=O)[O-]",
        "FADH(.)": "Cc1cc2[N]c3c([nH]c(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)OP(O)(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C"
    }
    
    for name, smi in test_smiles.items():
        flag, reason = is_flavin(smi)
        print(f"{name}: {flag} – {reason}")