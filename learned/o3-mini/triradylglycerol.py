"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol (functional parent glycerol CHEBI:17754, child: triglyceride CHEBI:17855, parent: glycerolipid CHEBI:35741)

A triradylglycerol is defined as a glycerol molecule in which each of the three possible positions (sn-1, sn-2, and sn-3) has one substituent group.
Each substituent group must be either acyl (typically bound as an ester), alkyl (bound as an ether) or alk-1-enyl (vinyl ether).
This script attempts to identify a glycerol core with three substituents and then attempts to classify each substituent.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    
    For this purpose we first search for a glycerol-like core defined as a chain
    of three carbons with the pattern CH2–CH–CH2. Then, for each of these carbons,
    we require that there is exactly one oxygen neighbor (not part of the core).
    The substituent attached at that oxygen is classified as:
      - acyl: if that substituent carbon shows a double bond to oxygen (i.e. a carbonyl),
      - alk-1-enyl: if that substituent carbon is sp2 hybridized,
      - alkyl: otherwise (assumed as a saturated alkyl chain).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a triradylglycerol, False otherwise
        str: Reason for classification or failure message
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # We search for a glycerol-like core.
    # The glycerol backbone in glycerol (and its substituted derivatives) is CH2–CH–CH2.
    # We use the SMARTS pattern for such a contiguous carbon chain.
    glycerol_smarts = "[CH2]-[CH]-[CH2]"
    glycerol_core = Chem.MolFromSmarts(glycerol_smarts)
    matches = mol.GetSubstructMatches(glycerol_core)
    if not matches:
        return False, "No glycerol backbone pattern found"

    # Iterate over candidate glycerol cores.
    for match in matches:
        # match is a tuple of atom indices for three contiguous carbons.
        # We assign positions: sn1, sn2 and sn3.
        core_atoms = [mol.GetAtomWithIdx(i) for i in match]
        substituents = {}  # dictionary to hold the substituent oxygen for each position
        valid_core = True  # flag to mark if we have a valid candidate glycerol core
        
        # Look for exactly one substituent oxygen (not part of the backbone) on each core carbon.
        for pos, atom in zip(['sn1','sn2','sn3'], core_atoms):
            oxy_neighbors = []
            for nb in atom.GetNeighbors():
                # Exclude atoms that belong to the core chain.
                if nb.GetAtomicNum() == 8 and nb.GetIdx() not in match:
                    oxy_neighbors.append(nb)
            if len(oxy_neighbors) != 1:
                valid_core = False
                # We note which position did not pass the expectation.
                break
            else:
                substituents[pos] = oxy_neighbors[0]
        if not valid_core:
            continue  # try next candidate glycerol core

        # For each substituent oxygen, identify the “other” neighbor which constitutes the head of the substituent group.
        substituent_types = {}
        for pos, oxy in substituents.items():
            # Find the substituent atom on the oxygen side (exclude the core carbon already connected)
            sub_atoms = [nb for nb in oxy.GetNeighbors() if nb.GetIdx() not in match]
            if not sub_atoms:
                valid_core = False
                break
            sub_atom = sub_atoms[0]  # assume first neighbor; in many cases the oxygen is attached to a single substituent group

            # Identify type
            # Check if the substituent group is acyl: 
            # Look for a double bond from sub_atom to an oxygen (carbonyl group) i.e. pattern O=C
            is_acyl = False
            for nb in sub_atom.GetNeighbors():
                # GetBondBetweenAtoms returns the bond object; bond type DOUBLE indicates double bond.
                bond = mol.GetBondBetweenAtoms(sub_atom.GetIdx(), nb.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nb.GetAtomicNum() == 8:
                    is_acyl = True
                    break
            if is_acyl:
                substituent_types[pos] = "acyl"
            else:
                # Check if the substituent carbon is sp2 hybridized, indicating an alk-1-enyl (vinyl ether) group.
                if sub_atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                    substituent_types[pos] = "alk-1-enyl"
                else:
                    # Otherwise, we classify it as alkyl.
                    substituent_types[pos] = "alkyl"
        if not valid_core or len(substituent_types) != 3:
            continue

        # If the candidate core has three substituents that can be classified as acyl, alkyl, or alk-1-enyl,
        # then we consider the molecule a triradylglycerol.
        reason = f"Found glycerol backbone with substituents: {substituent_types}"
        return True, reason

    return False, "No valid glycerol backbone with three proper substituents found"

# Example usage:
if __name__ == "__main__":
    # Try one of the provided SMILES strings
    smiles_example = "O(C(=O)CCCCCCCCCCCCCCCCC)C[C@H](OC(=O)CCCCCCC/C=C\\CCCC)COC(=O)CCCCCCCCCCCCCC"
    is_tri, msg = is_triradylglycerol(smiles_example)
    print(is_tri, msg)