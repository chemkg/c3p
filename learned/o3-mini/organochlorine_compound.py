"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: Organochlorine compound
Definition: An organochlorine compound is defined (here) as an organic molecule that 
contains at least one carbon–chlorine (C–Cl) bond that is considered chemically relevant.
Because (a) sometimes very large molecules may acquire a solitary Cl atom that is incidental 
and (b) some molecules have many Cl substituents, we use a two‐step heuristic:
  (1) first we look for any C–Cl bond;
  (2) then if found we check the overall chlorine-to‐carbon ratio and the molecular weight.
For small molecules (here, <500 Da) any C–Cl bond is taken as a positive hit.
For larger molecules, a single C–Cl bond is accepted only if the overall Cl/C ratio is at 
least 0.03 (i.e. about 3% chlorination); otherwise it is judged to be incidental.
If no direct C–Cl bond is found, we also allow for “highly chlorinated” molecules 
when the Cl/C ratio is very high (≥20%).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    
    The heuristic is as follows:
      1. Parse the molecule and count chlorine (atomic number 17) and carbon (atomic number 6) atoms.
      2. Look for at least one direct carbon–chlorine bond.
      3. Compute the chlorine-to–carbon ratio and the molecular weight.
         – For small molecules (MW < 500 Da): any clear C–Cl bond is taken as defining.
         – For larger molecules, a single C–Cl bond will only be accepted if the Cl/C ratio 
           is at least 0.03 (about 3%), to try to avoid incidental chlorination.
      4. In the (unlikely) event no direct C–Cl bond is found, we consider the molecule 
         “chlorinated” only if the ratio is very high (≥20%).
         
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
        
    # Look for a direct C–Cl bond.
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
    
    # Heuristic decision:
    if found_c_cl:
        # If there is only one chlorine atom:
        if len(cl_atoms) == 1:
            # For small molecules, assume that one C–Cl bond is significant.
            if mol_wt < 500:
                return True, "Molecule is small and contains a direct carbon–chlorine bond."
            else:
                # For larger molecules, require at least a minimal chlorine fraction.
                if cl_to_c_ratio >= 0.03:
                    return True, ("Molecule is large and the single C–Cl bond comes with a chlorine/carbon ratio "
                                  "of {:.2f} (≥0.03), indicating a likely organochlorine compound.").format(cl_to_c_ratio)
                else:
                    return False, ("Found one C–Cl bond but the chlorine content is very low (ratio: {:.2f}) in a large molecule; "
                                   "likely an incidental substitution.").format(cl_to_c_ratio)
        else:
            # If two or more Cl atoms (and at least one is directly bonded to C), accept it.
            return True, "Molecule contains multiple carbon–chlorine bonds."
    else:
        # If no explicit C–Cl bond is found, check if the chlorine fraction is extremely high.
        if cl_to_c_ratio >= 0.20:
            return True, ("No direct C–Cl bond was found, but the chlorine to carbon ratio is very high "
                          "({:.2f}), suggesting a highly chlorinated organic compound.").format(cl_to_c_ratio)
        else:
            return False, "No direct carbon–chlorine bonds found and chlorine content is too low."

# Example usage:
if __name__ == '__main__':
    # Test with a few examples:
    examples = [
        "OC(=O)\\C=C/C=C(/Cl)C(O)=O",  # 2-chloro-cis,cis-muconic acid (should be True)
        "COC1=CC=C(C=C1)C=NNC(=O)C2=CC=C(O2)COC3=CC=CC=C3Cl",  # a false-positive candidate
        "CC(=O)CCl"  # simple chloro-ketone
    ]
    for smi in examples:
        result, reason = is_organochlorine_compound(smi)
        print("SMILES:", smi)
        print("Classified as organochlorine?", result)
        print("Reason:", reason)
        print("----")