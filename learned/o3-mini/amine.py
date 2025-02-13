"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
A compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.
This implementation scans for nitrogen atoms that (a) are non‐aromatic and sp³–hybridized and then (b) have a total substituent
count (heavy neighbors + implicit hydrogens) consistent with a primary, secondary or tertiary amine. In the case of a positively charged
nitrogen (e.g. in a quaternary ammonium) a total count of 4 is allowed. Finally, if the nitrogen is attached to a carbon, we check whether that carbon 
is “carbonyl‐like” (i.e. double‐bonded to oxygen). Such N–C(=O) groups are generally considered amide‐like and are excluded EXCEPT when the carbonyl
is “thio‐modified” (i.e. the carbon is also bound to at least one sulfur).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_amine(smiles: str):
    """
    Determines if a molecule contains an amine functional group (derived from ammonia).
    
    The procedure is as follows:
      1. The SMILES is parsed to an RDKit molecule.
      2. For each nitrogen atom (atomic number 7) in the molecule:
          • It must be non‐aromatic and have sp³ hybridization.
          • Its total substituent count (explicit heavy atom neighbors + implicit hydrogens) should equal 3
            (as in a primary, secondary or tertiary amine) – or if the nitrogen bears a positive formal charge then we allow a count of 4.
          • We then check whether any neighbor (if it is a carbon) is “carbonyl–like” – that is, the carbon has at least one double
            bond to oxygen. However, if that carbon also is bonded (by a single bond) to a sulfur atom, then we do not exclude the candidate.
      3. If any nitrogen meets these criteria, we return True along with an explanatory reason.
    
    Args:
        smiles (str): A SMILES string for the molecule under investigation.
        
    Returns:
        (bool, str): A tuple (True, reason) if at least one amine group is found,
                     otherwise (False, reason).
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms looking for nitrogen candidates.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        
        # Exclude aromatic N atoms – many aromatic heterocycles or conjugated systems are not our intended amine
        if atom.GetIsAromatic():
            continue
        
        # Check that nitrogen is sp3-hybridized.
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        
        # Calculate the total substituent count: neighbors (heavy atoms) plus implicit hydrogens.
        tot_subs = len(atom.GetNeighbors()) + atom.GetTotalNumHs()
        
        # For a neutral amine (primary, secondary or tertiary) we expect tot_subs == 3.
        # If the nitrogen is protonated (formal charge > 0) then often tot_subs == 4.
        if (atom.GetFormalCharge() == 0 and tot_subs != 3) or (atom.GetFormalCharge() > 0 and tot_subs != 4):
            continue
        
        # Now check for "amide-like" attachment.
        # For each neighbor that is a carbon, look for a double bond to oxygen.
        reject_candidate = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                for bond in nbr.GetBonds():
                    # Check if the bond is a double bond.
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            # We have found a carbonyl-like bond.
                            # Now check if this carbon (nbr) is also bound to at least one sulfur (by a single bond).
                            has_sulfur = any(n.GetAtomicNum() == 16 for n in nbr.GetNeighbors() if n.GetIdx() != atom.GetIdx())
                            if not has_sulfur:
                                reject_candidate = True
                                break
                if reject_candidate:
                    break
        
        if reject_candidate:
            # Skip the candidate if it sits next to an unmodified carbonyl group.
            continue
        
        # If we reached here, we have found a nitrogen that meets our criteria.
        return True, "Found at least one amine functional group derived from ammonia"
    
    return False, "No amine functional group (derived from ammonia) identified"

# Uncomment below to do a test run with one of the provided examples.
# test_smiles = "CCCN(CCCC)CC(O)c1cc(Cl)cc2\\C(=C/c3ccc(Cl)cc3)c3cc(Cl)ccc3-c12"  # lumefantrine example
# print(is_amine(test_smiles))

# For further testing, one can loop over a list of SMILES strings.