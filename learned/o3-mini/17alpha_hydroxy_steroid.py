"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 17alpha-hydroxy steroid
Definition: The alpha-stereoisomer of 17-hydroxy steroid.
Heuristic:
  1. The molecule must have a steroid-like fused ring system. As a simple proxy,
     we require that the molecule has exactly one five-membered ring (ring D) and
     at least three six-membered rings (rings A, B, C). 
  2. Within the five-membered ring, we search for a carbon atom that is substituted by a hydroxyl group (-OH).
  3. We then require that this carbon is marked as chiral and that the CIP stereochemistry,
     as assigned by RDKit, is “S” (which here we assume corresponds to the α configuration).
Note: This is a heuristic approach. Stereochemistry in steroids is complex and a more
rigorous approach may require a specific mapping of the steroid numbering.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule belongs to the 17alpha-hydroxy steroid class based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 17alpha-hydroxy steroid, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is computed (including CIP codes)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Get ring information for the molecule.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings detected – not a steroid"
    
    # Identify five-membered rings; steroid D-ring is five-membered.
    five_membered = [ring for ring in rings if len(ring) == 5]
    if len(five_membered) != 1:
        return False, f"Expected 1 five-membered ring (steroid D-ring), but found {len(five_membered)}"
    
    # For extra confidence, check that there are at least three six-membered rings.
    six_membered = [ring for ring in rings if len(ring) == 6]
    if len(six_membered) < 3:
        return False, f"Expected at least 3 six-membered rings (steroid A, B, C), but found {len(six_membered)}"
    
    # Define a SMARTS to find hydroxyl groups: –OH (oxygen with an attached hydrogen)
    oh_smarts = Chem.MolFromSmarts("[OX2H]")
    
    # For the one 5-membered ring, consider each atom.
    d_ring = five_membered[0]
    candidate_found = False
    candidate_reason = ""
    
    for atom_idx in d_ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        # We require a carbon (ring D is made of carbons in a steroid nucleus)
        if atom.GetAtomicNum() != 6:
            continue
        # Check if the atom has an –OH group attached
        neighbors = atom.GetNeighbors()
        has_oh = False
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 8:
                # Use the SMARTS match to check for hydroxyl oxygen
                if nbr.HasSubstructMatch(oh_smarts) or nbr.GetTotalNumHs() > 0:
                    has_oh = True
                    break
        if not has_oh:
            continue

        # Check if this carbon is marked as chiral.
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            candidate_reason = "Hydroxyl-bearing carbon in D-ring is not stereochemically defined"
            continue
        
        # Check if CIP stereochemistry is assigned.
        if not atom.HasProp('_CIPCode'):
            candidate_reason = "No CIP stereochemistry assigned for hydroxyl-bearing carbon in D-ring"
            continue
        cip = atom.GetProp('_CIPCode')
        # In this heuristic, we assume that "S" corresponds to the alpha configuration.
        if cip == "S":
            candidate_found = True
            candidate_reason = "Found steroid nucleus with D-ring hydroxyl carbon having CIP 'S' (α)"
            break
        else:
            candidate_reason = f"D-ring hydroxyl carbon has CIP '{cip}' which does not match expected 'S' for alpha"
    
    if candidate_found:
        return True, candidate_reason
    else:
        # If no candidate in the D-ring satisfies the requirement, return False.
        return False, f"Steroid nucleus identified but no 17α–OH group detected. {candidate_reason}"

# Example usage (you can remove or comment these lines when using this module in production):
if __name__ == "__main__":
    # Example SMILES for 17α-hydroxypregnenolone (one of the listed examples)
    test_smiles = "CC(=O)[C@@]1(O)CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_17alpha_hydroxy_steroid(test_smiles)
    print(result, reason)