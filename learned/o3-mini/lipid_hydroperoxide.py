"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide
Defined as any lipid carrying one or more hydroperoxy (-OOH) substituents.
A lipid hydroperoxide in this context must contain at least one hydroperoxy group 
and possess lipid-like characteristics (e.g. a long aliphatic chain indicated by a minimum
number of carbon atoms and an appropriate molecular weight).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    
    A lipid hydroperoxide is defined as a lipid (here, judged by a long carbon chain and a certain
    molecular weight) that carries one or more hydroperoxy substituents (-OOH).
    
    We use the following heuristic tests:
      1. The molecule must parse correctly.
      2. It must have at least one hydroperoxy group. The SMARTS "[OX2]-[OX2H]" is used to detect
         an -O-OH fragment (the first oxygen atom connected to a second oxygen that bears one hydrogen).
      3. The molecule should be lipid-like. Here we require a minimum of 12 carbon atoms.
      4. The molecular weight should be at least 250 Da, a rough lower bound for typical lipids.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lipid hydroperoxide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that the hydroperoxy O with its H can be detected.
    mol = Chem.AddHs(mol)

    # Identify the hydroperoxy group (-OOH).
    # The SMARTS "[OX2]-[OX2H]" attempts to capture an oxygen (with two substituents) bonded to
    # another oxygen that carries one hydrogen.
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2]-[OX2H]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy (-OOH) substituent found"

    # Count the number of carbons in the molecule.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Insufficient number of carbon atoms ({carbon_count}) for a lipid-like structure"

    # Calculate the molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 250:
        return False, f"Molecular weight too low ({mw:.1f} Da) for a typical lipid"

    return True, "Molecule contains a hydroperoxy substituent and displays lipid-like characteristics"