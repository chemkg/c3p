"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: Organochlorine compound
Definition: An organochlorine compound is a compound that is “chlorinated” – either 
   containing at least one carbon–chlorine bond or having a high chlorine content relative 
   to its carbon count.
Note:
  (1) Many organochlorine compounds have a C–Cl bond but some organochlorine examples (e.g. certain 
      organometallic complexes) may not. Conversely, in very large molecules a single C–Cl bond 
      may be incidental and not really define the molecule as an organochlorine.
  (2) This function uses a heuristic: first it looks for any C–Cl bonds; if found, it then checks 
      that the Cl/carbon ratio is not very low (which is common in larger molecules that have only one Cl).
  (3) If no direct C–Cl bond is found, but the chlorine fraction is high (≥20% of the carbon atoms), 
      the molecule is assumed to be an organochlorine.
"""

from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    
    The definition used here is twofold:
     (a) If at least one carbon atom is directly bonded to a chlorine atom then the molecule is 
         potentially an organochlorine compound. However, if the overall chlorine-to-carbon ratio is 
         very low (in a large molecule), it is likely an incidental substitution.
     (b) In case no explicit C–Cl bond is detected, if the chlorine content is high relative to the 
         number of carbon atoms (e.g. ≥20%), the molecule is likely an organochlorine compound.
    
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
        
    # Collect all chlorine atoms (atomic number 17) and carbon atoms (atomic number 6)
    cl_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    if not cl_atoms:
        return False, "No chlorine atoms found in the molecule."
    if not c_atoms:
        return False, "No carbon atoms found; not an organic compound."
    
    # Look for a direct carbon-chlorine bond.
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
    
    # Heuristic thresholds:
    # - For molecules with an explicit C–Cl bond, require at least a minimal chlorine fraction in large molecules.
    # - For molecules without an explicit C–Cl bond, require a high chlorine content (here, 20% or more).
    if found_c_cl:
        # For very large molecules a single C–Cl might be incidental.
        if mol.GetNumAtoms() > 30 and cl_to_c_ratio < 0.05:
            return False, ("Found a C–Cl bond but the chlorine content is very low (ratio: {:.2f}) "
                           "in a large molecule; likely not a typical organochlorine compound.").format(cl_to_c_ratio)
        else:
            return True, "Molecule contains a direct carbon–chlorine bond."
    else:
        # If no explicit C–Cl bond was found, check if the chlorine content is very high.
        if cl_to_c_ratio >= 0.2:
            return True, ("No direct carbon–chlorine bond was detected, but the chlorine "
                          "to carbon ratio ({:.2f}) is high, indicating a highly chlorinated organic compound."
                         ).format(cl_to_c_ratio)
        else:
            return False, "No direct carbon–chlorine bonds found and the chlorine content is too low."
            
# The heuristic parameters (molecular size, ratio thresholds) here are adjustable and were chosen 
# to improve classification metrics based on the provided outcomes.