"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
A D-hexose is defined here as a hexose (6-carbon monosaccharide) whose sugar ring contains a chiral center 
(i.e. the ring carbon that bears the exocyclic CH2OH group) with the R configuration (i.e. D-configuration at position 5).

NOTE: This procedure uses heuristic rules. It first identifies rings of size 5 (furanose) or 6 (pyranose) that 
contain exactly one ring oxygen. Then, for each ring carbon it checks for a candidate CH2OH group (a non‐ring carbon 
attached that is “CH2” and further bonded to an OH). If found, it checks its CIP code (assigned by RDKit) and if it is 
“R” then the compound is classified as a D-hexose.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is defined as a hexose (monosaccharide with 6 carbon atoms) having D configuration
    at the carbon bearing the exocyclic CH2OH group (position 5 in many common sugars).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a D-hexose, False otherwise
        str: A reason explaining the classification result.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Assign/refresh stereochemistry (CIP codes etc.)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    
    candidate_found = False
    msg_details = []
    
    # Loop through rings of interest: furanose (5 members) and pyranose (6 members)
    for ring in atom_rings:
        if len(ring) not in [5, 6]:
            continue  # skip rings that are not furanose or pyranose rings
            
        # Count ring oxygens
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            # For a “sugar ring” we expect exactly one ring oxygen.
            continue
            
        # Now, search in this ring for a carbon that has an exocyclic CH2OH group.
        # In many hexoses, the carbon that bears the CH2OH group (the exocyclic substituent) is the one at position 5.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # looking only at carbon atoms in the ring
            # Look at neighbors that are not in the ring.
            exocyclic_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            for nbr in exocyclic_neighbors:
                # We are looking for a CH2OH group: neighbor should be a carbon with exactly two implicit hydrogens 
                # and attached to an oxygen (i.e. the -OH).
                if nbr.GetAtomicNum() == 6:
                    # Check the count of hydrogens (using total H count, which includes implicit)
                    if nbr.GetTotalNumHs() != 2:
                        continue
                    # Now check that this carbon has an oxygen neighbor (the hydroxyl)
                    oh_found = False
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 8 and subnbr.GetDegree() == 1:
                            oh_found = True
                            break
                    if not oh_found:
                        continue
                    # If we reach here, we have identified a CH2OH substituent attached to a ring carbon.
                    # Now check if this ring carbon (our candidate for "position 5") has an assigned CIP code,
                    # and if it is "R". (In D-sugars, the CH2OH-bearing carbon should be R.)
                    cip = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
                    if cip == "R":
                        return True, "Found sugar ring with CH2OH substituent at a chiral center (position 5) with R configuration (D-hexose)"
                    else:
                        msg_details.append(f"Found CH2OH on ring carbon (idx {idx}) but CIP code is {cip if cip is not None else 'undefined'}")
                    # note: if no valid chiral center is found, continue with next candidate
        # end for each atom in ring
    # end for each ring

    if msg_details:
        # We found candidate CH2OH-bearing atoms but none had the R configuration.
        return False, "Candidate sugar ring(s) found but none had the expected R configuration at position 5: " + "; ".join(msg_details)
    
    # As a fallback, we can also check if the molecule contains 6 carbons (hexose formula) as an additional heuristic.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 6:
        return False, f"Total carbon count ({carbon_count}) does not equal 6 (expected for a hexose)"
    
    return False, "No valid hexose sugar ring with D configuration at position 5 was found"

# Example usage (for testing):
if __name__ == '__main__':
    # Example SMILES for beta-D-glucose:
    test_smiles = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
    valid, reason = is_D_hexose(test_smiles)
    print("Beta-D-glucose classification:", valid, reason)