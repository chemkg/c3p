"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
A compound formally derived from ammonia by replacing one, two or three hydrogen atoms 
by hydrocarbyl groups.

This implementation inspects each nitrogen atom in the molecule and:
  1. Accepts both sp2 and sp3 nitrogen if it is not embedded in an aromatic ring.
  2. Checks that the total substituent count (explicit heavy neighbors + implicit hydrogens)
     equals 3 for neutral amines (as in primary, secondary or tertiary amines) or 4 for 
     positively charged (e.g. quaternary) ammoniums.
  3. Rejects a candidate if any neighboring carbon is “carbonyl‐like” (i.e. has a double bond 
     to oxygen) and that carbon is not “thio‐modified” (i.e. lacking a single-bonded sulfur neighbor).
     
This design aims to reduce both the false negatives (by accepting sp2 amines not embedded
in a ring) and the false positives (by excluding nitrogens adjacent to “amide‐like” carbonyls).
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_amine(smiles: str):
    """
    Determines if a molecule contains an amine functional group (derived from ammonia).
    
    The procedure is as follows:
      1. Parse the SMILES to an RDKit molecule.
      2. For each nitrogen (atomic number 7):
           - Skip if the nitrogen is part of an aromatic ring (e.g. pyridine) 
             because such ring‐members are generally not considered ammonia derivatives.
           - Compute the total substituent count as the number of heavy atom neighbors plus implicit Hs.
           - For a neutral amine, require a count of 3 (as in RNH2, R2NH, or R3N).
             For a positively charged nitrogen, require a count of 4.
           - For each neighbor carbon, if that carbon bears a double bond to oxygen,
             then check if that same carbon also is single‐bonded to at least one sulfur.
             If not, then the candidate is likely an amide and is rejected.
      3. Return True (plus a reason) as soon as one candidate passes these criteria.
         If none are found, return False.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        (bool, str): A tuple: True and a reason if at least one amine group is found;
                     otherwise False and an explanatory reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop over atoms, searching for candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue

        # Exclude nitrogen atoms that are members of an aromatic ring.
        # (This avoids, for example, pyridine-like nitrogens.)
        if atom.IsInRing() and atom.GetIsAromatic():
            continue

        # Calculate the total substituent count:
        # (number of heavy atom neighbors + implicit hydrogens).
        tot_subs = len(atom.GetNeighbors()) + atom.GetTotalNumHs()
        fc = atom.GetFormalCharge()
        # For a neutral amine, we expect three substituents
        if fc == 0 and tot_subs != 3:
            continue
        # For a positively charged nitrogen (e.g. quaternary ammonium), permit a total of 4.
        if fc > 0 and tot_subs != 4:
            continue

        # Check: if the nitrogen is directly attached to a carbonyl-like carbon, then exclude.
        reject_candidate = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:  # oxygen found (carbonyl)
                            # Check if this carbon (nbr) is also attached to at least one sulfur (by a single bond).
                            has_sulfur = any(
                                (b.GetBondType() == rdchem.BondType.SINGLE and 
                                 b.GetOtherAtom(nbr).GetAtomicNum() == 16)
                                for b in nbr.GetBonds()
                            )
                            if not has_sulfur:
                                reject_candidate = True
                                break
                if reject_candidate:
                    break
        
        if reject_candidate:
            continue

        # If we reached here, we have found an appropriate amine candidate.
        return True, "Found at least one amine functional group derived from ammonia"
    
    return False, "No amine functional group (derived from ammonia) identified"

# Example (uncomment to test with one of the provided SMILES):
# test_smiles = "CCCN(CCCC)CC(O)c1cc(Cl)cc2\\C(=C/c3ccc(Cl)cc3)c3cc(Cl)ccc3-c12"  # lumefantrine example
# print(is_amine(test_smiles))