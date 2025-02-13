"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: Beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
This implementation attempts to detect a beta-D-glucoside moiety by first identifying a
pyranose (6-membered sugar) ring that:
  • Contains exactly 1 ring oxygen and 5 carbons.
  • Has at least one ring carbon that bears an exocyclic oxygen substituent.
  • That exocyclic oxygen is attached to a CH2 group (i.e. –CH2OH), enforcing the presence
    of the characteristic exocyclic CH2OH of D-glucosides.
  • The candidate anomeric carbon is stereochemically defined (a proxy for having beta configuration).
If such a ring is found, the molecule is classified as a beta-D-glucoside.
Note: Stereochemistry in sugar chemistry is subtle and SMILES representations vary widely.
Thus this algorithm is an approximation that aims to improve on the previous pure SMARTS approach.
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    Rather than rely solely on strict SMARTS patterns, this function uses an algorithmic
    approach:
      1. Parse the molecule.
      2. Identify 6-membered rings (pyranoses) that contain exactly 1 oxygen and 5 carbons.
      3. For each such ring, look for a candidate “anomeric” carbon – a carbon in the ring
         that is exocyclically bound to an oxygen. In a beta-D-glucoside this oxygen should
         lead to a CH2 group (i.e. be part of a CH2OH group).
      4. Check that the candidate carbon has its chirality specified (i.e. not CHI_UNSPECIFIED).
         This is a proxy for having the beta configuration.
         
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a beta-D-glucoside fragment is detected, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that ring information is updated
    mol.UpdatePropertyCache()
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    for ring in rings:
        # We are only interested in 6-membered rings (pyranoses).
        if len(ring) != 6:
            continue
        
        # Count number of ring oxygen and carbons.
        oxygens_in_ring = 0
        carbons_in_ring = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygens_in_ring += 1
            elif atom.GetAtomicNum() == 6:
                carbons_in_ring += 1
        # For a glucopyranose, we expect exactly 1 oxygen and 5 carbons.
        if oxygens_in_ring != 1 or carbons_in_ring != 5:
            continue
        
        # Now, search for an "anomeric" carbon candidate within the ring.
        # The anomeric carbon is usually a ring carbon that has an exocyclic oxygen.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Skip if not carbon.
            if atom.GetAtomicNum() != 6:
                continue
            
            # Look for exocyclic oxygen neighbors (neighbors not in the ring).
            exo_oxys = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum()==8 and nbr.GetIdx() not in ring]
            if not exo_oxys:
                continue
            
            # To enforce that the exocyclic oxygen is part of a CH2OH group,
            # check if at least one exocyclic oxygen is bound to a carbon that has two hydrogens.
            ch2oh_found = False
            for oxy in exo_oxys:
                for nbr in oxy.GetNeighbors():
                    # Exclude the link back to our candidate anomeric carbon.
                    if nbr.GetIdx() == atom.GetIdx():
                        continue
                    # Check if this neighbor is carbon and has two attached hydrogens.
                    if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 2:
                        ch2oh_found = True
                        break
                if ch2oh_found:
                    break
            if not ch2oh_found:
                continue
            
            # Check that the candidate carbon has defined chirality.
            # If the chiral tag is CHI_UNSPECIFIED, we cannot tell beta configuration.
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                continue
            
            # If we have reached here, we have a candidate pyranose ring with an anomeric carbon
            # that is both exocyclicly linked (with a CH2OH group) and stereochemically defined.
            return True, "Beta-D-glucoside fragment detected in a pyranose ring."
    
    # If no candidate ring was found
    return False, "No beta-D-glucoside fragment found in the molecule."

# Example usage for testing:
if __name__ == "__main__":
    # Test on one example from the provided list:
    test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
    result, reason = is_beta_D_glucoside(test_smiles)
    print("Result:", result)
    print("Reason:", reason)