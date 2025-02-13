"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: Organochlorine compound
Definition: An organochlorine compound is defined (here) as an organic molecule 
that contains at least one carbon–chlorine (C–Cl) bond that is deemed chemically relevant.
Because very large molecules sometimes have a solitary Cl that is incidental, we use the following heuristic:
  1. Parse the molecule and count chlorine (atomic number 17) and carbon (atomic number 6) atoms.
  2. Look for at least one explicit C–Cl bond.
  3. Compute the chlorine-to–carbon ratio and the molecular weight.
     – For small molecules (MW < 500 Da): if a single Cl is found, require the Cl/C ratio to be at least 0.05.
     – For larger molecules (MW ≥ 500 Da): a single C–Cl bond is accepted only if the Cl/C ratio is at least 0.03.
  4. In any case, if two or more chlorine atoms are directly bonded, accept the molecule.
  5. Also, if no direct C–Cl bond is found we allow for “highly chlorinated” molecules if the Cl/C ratio is very high (≥ 0.20). 

This updated heuristic aims to reduce false positives (molecules poorly classified as organochlorine) 
while still capturing the intended molecules.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    
    The heuristic is as follows:
      1. Parse the molecule and count chlorine (atomic number 17) and carbon (atomic number 6) atoms.
      2. Check if at least one direct C–Cl bond exists.
      3. Compute the chlorine-to–carbon ratio and the molecular weight.
         • For small molecules (MW < 500 Da): 
              - If there is only one Cl, require the ratio to be at least 0.05.
         • For large molecules (MW >= 500 Da): 
              - A single C–Cl bond is accepted only if the Cl/C ratio is at least 0.03.
      4. If two or more Cl atoms are present (with at least one direct C–Cl bond), accept the molecule.
      5. If no direct C–Cl bond is detected, then only accept if the Cl/C ratio is extremely high (≥ 0.20).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organochlorine compound, False otherwise.
        str: Explanation of the result.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Count chlorine and carbon atoms.
    cl_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    c_atoms  = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    if not cl_atoms:
        return False, "No chlorine atoms found in the molecule."
    if not c_atoms:
        return False, "No carbon atoms found in the molecule; not an organic compound."
        
    # Look for a direct C–Cl bond by checking neighbors of chlorine atoms.
    found_c_cl = False
    for cl in cl_atoms:
        for neighbor in cl.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                found_c_cl = True
                break
        if found_c_cl:
            break

    # Compute the ratio of chlorine atoms to carbon atoms.
    cl_to_c_ratio = len(cl_atoms) / len(c_atoms)
    
    # Compute the molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Use the updated heuristic:
    if found_c_cl:
        # If two or more chlorine atoms are directly bonded, accept.
        if len(cl_atoms) >= 2:
            return True, "Molecule contains multiple carbon–chlorine bonds."
        else:
            # Exactly one chlorine atom is directly bonded.
            if mol_wt < 500:
                # For small molecules, require a minimum Cl/C ratio (here set to 0.05).
                if cl_to_c_ratio >= 0.05:
                    return True, "Molecule is small and has a significant direct C–Cl bond (Cl/C ratio >= 0.05)."
                else:
                    return False, ("Molecule is small and contains one C–Cl bond, "
                                   "but the chlorine-to-carbon ratio (%.2f) is too low." % cl_to_c_ratio)
            else:
                # For larger molecules, require a Cl/C ratio of at least 0.03.
                if cl_to_c_ratio >= 0.03:
                    return True, ("Molecule is large and the single C–Cl bond comes with a chlorine/carbon ratio "
                                  "of %.2f (>= 0.03), indicating a likely organochlorine compound." % cl_to_c_ratio)
                else:
                    return False, ("Found one C–Cl bond but the chlorine content is very low (ratio: %.2f) in a large molecule; "
                                   "likely an incidental substitution." % cl_to_c_ratio)
    else:
        # No explicit C–Cl bond detected. Accept only if the entire molecule is extremely chlorinated.
        if cl_to_c_ratio >= 0.20:
            return True, ("No direct C–Cl bond was found, but the chlorine to carbon ratio is very high "
                          "(%.2f), suggesting a highly chlorinated organic compound." % cl_to_c_ratio)
        else:
            return False, "No direct carbon–chlorine bonds found and chlorine content is too low."

# Example usage:
if __name__ == '__main__':
    # List of test examples (both those expected to be positive and some known false positives)
    examples = [
        "OC(=O)\\C=C/C=C(/Cl)C(O)=O",  # 2-chloro-cis,cis-muconic acid: expected positive (high ratio)
        "C1C[C@@H]2C(C(=C([C@H]1C2)O)C(C3=CC=C(C=C3Cl)S(C)(=O)=O)=O)=O",  # benzobicyclon hydrolysate: expected positive
        "C1=CC=C(C(=C1)CN2C=CSC2=NC(=O)CCl)Cl",  # 2-chloro-N-[3-[(2-chlorophenyl)methyl]-2-thiazolylidene]acetamide: expected positive
        "COC1=CC=C(C=C1)C=NNC(=O)C2=CC=C(O2)COC3=CC=CC=C3Cl",  # false-positive candidate from previous attempt
        "CC(=O)CCl"  # simple chloro ketone
    ]
    for smi in examples:
        result, reason = is_organochlorine_compound(smi)
        print("SMILES:", smi)
        print("Classified as organochlorine?", result)
        print("Reason:", reason)
        print("----")