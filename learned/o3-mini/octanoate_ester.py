"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: Octanoate ester
Definition: Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).

Improvement rationale:
  • Instead of simply verifying the presence of an ester, we inspect each ester group found.
  • For each ester, we identify the acyl (acid) side and “walk” along its carbon skeleton
    following only single bonds between sp3 carbons, ensuring that the chain is linear.
  • We then check that exactly eight carbons are present (the carbonyl carbon plus 7 more)
    and that the terminal carbon is a methyl group.
  • Only then do we consider the ester to be derived from octanoic acid.
  
Note: This method may not be perfect for every edge case but is a step toward reducing false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def _get_linear_acyl_chain(mol: rdchem.Mol, carbonyl, acyl_neighbor):
    """
    Traverses from the carbonyl carbon through the acyl chain
    following single bonds between sp3 carbons.
    
    Args:
        mol (rdchem.Mol): The molecule.
        carbonyl (Atom): The carbonyl carbon atom (starts the chain).
        acyl_neighbor (Atom): The neighbor (sp3 carbon) attached to the carbonyl that begins the chain.
        
    Returns:
        list: A list of atoms in order from the carbonyl to the end of the chain,
              if the chain is linear. Returns None if a branch is encountered.
    """
    # Start chain with the carbonyl carbon and the neighbor.
    chain = [carbonyl, acyl_neighbor]
    prev_atom = carbonyl
    current = acyl_neighbor
    
    while True:
        # Look for next carbon neighbors (exclude the atom we came from)
        next_atoms = []
        for nbr in current.GetNeighbors():
            if nbr.GetIdx() == prev_atom.GetIdx():
                continue
            # Only follow carbon atoms with single bonds.
            if nbr.GetAtomicNum() == 6:
                bond = mol.GetBondBetweenAtoms(current.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == rdchem.BondType.SINGLE:
                    next_atoms.append(nbr)
        if len(next_atoms) == 0:
            # End of chain reached.
            break
        elif len(next_atoms) == 1:
            # Continue linearly.
            prev_atom = current
            current = next_atoms[0]
            chain.append(current)
        else:
            # Branching encountered, so not a linear chain.
            return None
    return chain

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is defined as a fatty acid ester whose acid part is derived from octanoic acid.
    
    Octanoic acid has 8 carbons in total:
       CH3(CH2)6C(=O)OH
    In an ester the acidic –OH is replaced, so we want to see a chain of exactly 8 carbons
    (with the carbonyl counting as one) that is linear and with the terminal carbon being CH3.
    
    The algorithm:
      1. Parses the SMILES string.
      2. Finds ester groups using the SMARTS pattern "[#6X3](=O)[OX2H0]".
      3. For each match, identifies the acyl neighbor of the carbonyl
         (i.e. the carbon attached to the carbonyl, not the ester oxygen).
      4. Traverses the acyl chain in a strictly linear fashion.
      5. Checks that the acyl chain has exactly 8 carbons and ends in a methyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True and a reason if the molecule contains an octanoate ester moiety,
                     False and a reason otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for an ester group: a carbonyl carbon (sp2 carbon, atomic number 6, 3 connected atoms)
    # with a double bond to oxygen and a single bond to an oxygen (not an -OH)
    ester_smarts = "[#6X3](=O)[OX2H0]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return None, None
    
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "No ester group found in molecule"
    
    # Check each ester match.
    for match in matches:
        # In the SMARTS, match[0] is the carbonyl carbon and match[1] is the ester oxygen.
        carbonyl = mol.GetAtomWithIdx(match[0])
        ester_oxygen = mol.GetAtomWithIdx(match[1])
        
        # Identify the other neighbor of the carbonyl: the acyl chain side.
        acyl_neighbors = []
        for nbr in carbonyl.GetNeighbors():
            # Exclude the ester oxygen as well as the carbonyl's double-bonded oxygen.
            if nbr.GetIdx() == ester_oxygen.GetIdx():
                continue
            # The acyl chain should be a carbon.
            if nbr.GetAtomicNum() == 6:
                acyl_neighbors.append(nbr)
        if not acyl_neighbors:
            continue  # no acyl branch found in this ester
        # We check each possible acyl neighbor.
        for acyl in acyl_neighbors:
            # Get the full linear chain from the carbonyl outward.
            chain = _get_linear_acyl_chain(mol, carbonyl, acyl)
            if chain is None:
                continue  # chain was branched or not linear
            # We expect exactly 8 carbons in the acyl chain for octanoic acid.
            if len(chain) != 8:
                continue
            # Additionally, check the terminal carbon is a methyl group:
            # it should have exactly one carbon neighbor (the previous atom in the chain)
            terminal = chain[-1]
            carbon_neighbors = [nbr for nbr in terminal.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                continue
            # If all conditions are met, we classify this as an octanoate ester.
            return True, "Molecule contains an octanoate ester moiety"
        
    return False, "Molecule does not contain an octanoate ester moiety"

# Example usage and testing:
if __name__ == "__main__":
    test_smiles_list = [
        "CCCCCCCC(=O)OCC",  # ethyl octanoate - expect True
        "CCCCCCCC(=O)OC",   # methyl octanoate - expect True
        "CCCCCCCC(=O)OC[C@H](O)CO", # 1-octanoyl-sn-glycerol - expect True
        "CC(=O)OCCCCCCCC(=O)O",  # one octanoate, one acetyl ester - expect True only if an octanoate is present
        # A false positive example from prior results may be a phosphocholine ester (if its acyl chain is octanoyl)
        "P(OC[C@@H](CO)OC(CCCCCCC)=O)(=O)(OCC[N+](C)(C)C)[O-]",  # 2-octanoyl-sn-glycero-3-phosphocholine; expecting False
    ]
    for smi in test_smiles_list:
        result, reason = is_octanoate_ester(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")