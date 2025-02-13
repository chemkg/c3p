"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol
Defined as a glycerol compound where each of the three hydroxyl positions (sn-1, sn-2, sn-3)
is substituted with one group that is either acyl (ester, if its immediate carbon bears a carbonyl),
alkyl (saturated ether) or alk-1-enyl (vinyl ether).
The feature parent is glycerol (CHEBI:17754) and one child is triglyceride (CHEBI:17855)
with parent glycerolipid (CHEBI:35741).

Improvements:
 - Instead of simply finding an oxygen on each core carbon, we now require that each oxygen 
   has exactly two neighbors (only the glycerol carbon and one substituent carbon). 
 - We also check that the substituent carbon is actually a carbon (atomic number 6), not, for example, phosphorus.
 - The substituent is classified as acyl if its carbon shows a carbonyl bond,
   as alk-1-enyl if that carbon is sp2-hybridized, and otherwise as alkyl.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    
    This function first searches for a glycerol-like core defined as a chain of three carbons with the pattern CH2–CH–CH2.
    For each carbon in the candidate core, it requires exactly one oxygen neighbor not in the core.
    In addition, it verifies that:
      - Each substituent oxygen has exactly two neighbors (one being the glycerol carbon).
      - The oxygen’s other neighbor (the start of the substituent) is a carbon (atomic number 6) 
        and is not part of any interfering groups (such as phosphate).
      - The substituent is classified as follows:
            • "acyl" if that carbon is bound via a double bond to an oxygen (carbonyl group),
            • "alk-1-enyl" if that carbon is sp2 hybridized,
            • "alkyl" otherwise.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a triradylglycerol, False otherwise.
        str: Reason for classification or failure message.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a glycerol-like backbone: CH2-CH-CH2.
    glycerol_smarts = "[CH2]-[CH]-[CH2]"
    glycerol_core = Chem.MolFromSmarts(glycerol_smarts)
    matches = mol.GetSubstructMatches(glycerol_core)
    if not matches:
        return False, "No glycerol backbone pattern found"

    # Iterate over candidate glycerol backbones
    for match in matches:
        # Use the matched atom indices as a candidate glycerol core.
        core_indices = set(match)
        core_atoms = [mol.GetAtomWithIdx(i) for i in match]
        substituents = {}  # To hold the oxygen atom attached to each sn-position
        valid_core = True
        
        # For each carbon in the glycerol core, ensure exactly one oxygen neighbor (not part of the core).
        for pos, atom in zip(['sn1', 'sn2', 'sn3'], core_atoms):
            oxy_neighbors = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 8 and nb.GetIdx() not in core_indices]
            if len(oxy_neighbors) != 1:
                valid_core = False
                break
            # Check that the oxygen atom is only involved in two bonds (glycerol C and its substituent)
            oxy = oxy_neighbors[0]
            if oxy.GetDegree() != 2:
                valid_core = False
                break
            substituents[pos] = oxy

        if not valid_core:
            continue  # Try next candidate glycerol core
        
        # Classify each substituent from the oxygen atom.
        substituent_types = {}
        for pos, oxy in substituents.items():
            # Get the substituent atom: the neighbor of oxygen that is not the glycerol carbon.
            sub_neighbors = [nb for nb in oxy.GetNeighbors() if nb.GetIdx() not in core_indices]
            if len(sub_neighbors) != 1:
                valid_core = False
                break            
            sub_atom = sub_neighbors[0]
            # Ensure the substituent atom is a carbon; if not, this candidate is not a simple triradylglycerol.
            if sub_atom.GetAtomicNum() != 6:
                valid_core = False
                break
            # Also, check that this substituent oxygen is not part of a phosphate linkage:
            # (i.e. its other neighbor should not be phosphorus)
            for nb in oxy.GetNeighbors():
                if nb.GetAtomicNum() == 15:
                    valid_core = False
                    break
            if not valid_core:
                break
            
            # Determine the substituent type:
            # Check for acyl: if the substituent carbon (sub_atom) is bonded via a double bond to an oxygen (carbonyl)
            is_acyl = False
            for nb in sub_atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(sub_atom.GetIdx(), nb.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nb.GetAtomicNum() == 8:
                    is_acyl = True
                    break
            if is_acyl:
                substituent_types[pos] = "acyl"
            else:
                # Check for alk-1-enyl (vinyl ether): if the substituent carbon is sp2 hybridized
                if sub_atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                    substituent_types[pos] = "alk-1-enyl"
                else:
                    substituent_types[pos] = "alkyl"
        if not valid_core or len(substituent_types) != 3:
            continue
        
        reason = f"Found glycerol backbone with substituents: {substituent_types}"
        return True, reason

    return False, "No valid glycerol backbone with three proper substituents found"

# Example usage:
if __name__ == "__main__":
    # Example: one of the provided SMILES strings.
    smiles_example = "O(C(=O)CCCCCCCCCCCCCCCCC)C[C@H](OC(=O)CCCCCCC/C=C\\CCCC)COC(=O)CCCCCCCCCCCCCC"
    result, message = is_triradylglycerol(smiles_example)
    print(result, message)