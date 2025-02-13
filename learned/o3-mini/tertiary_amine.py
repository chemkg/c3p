"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: A tertiary amine
A tertiary amine is defined as a compound formally derived from ammonia by replacing 
three hydrogen atoms by hydrocarbyl groups. That is, the molecule must contain at least one
nitrogen atom (atomic number 7) that is neutral, has no hydrogen attached and is bonded
to exactly three carbon-containing groups.
"""

from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule contains a tertiary amine group based on its SMILES string.
    A tertiary amine is defined here as a nitrogen (atomic number 7) that is formally neutral,
    has no directly attached hydrogens and is connected to three other atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule contains a tertiary amine group, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over atoms to find a candidate tertiary amine.
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen (atomic number 7)
        if atom.GetAtomicNum() != 7:
            continue
        
        # For a tertiary amine, the nitrogen should be neutral.
        if atom.GetFormalCharge() != 0:
            continue
        
        # The nitrogen should have no directly attached hydrogen atoms.
        # GetTotalNumHs() returns both implicit and explicit hydrogens.
        if atom.GetTotalNumHs() != 0:
            continue
            
        # The nitrogen should have three neighbors (i.e. three substituents).
        if atom.GetDegree() != 3:
            continue
        
        # Optional: check that each substituent contains at least one carbon.
        # This step helps assure that the substituents are "hydrocarbyl" groups.
        all_substituents_are_hydrocarbyl = True
        for neighbor in atom.GetNeighbors():
            # If a neighbor is not carbon, we allow it only if it is not a hydrogen.
            # But to be strict to the definition, we expect a hydrocarbyl group to contain at least one carbon.
            if neighbor.GetAtomicNum() != 6:
                all_substituents_are_hydrocarbyl = False
                break
        
        if not all_substituents_are_hydrocarbyl:
            continue
        
        # Found a nitrogen atom matching our tertiary amine criteria.
        return True, f"Found tertiary amine center at atom index {atom.GetIdx()}"
        
    # If loop finishes without a match, no tertiary amine center has been found.
    return False, "No tertiary amine center found"

# Below are some simple tests (these lines can be removed if using this function as a module):
if __name__ == "__main__":
    examples = {
        "benzphetamine": "C[C@@H](Cc1ccccc1)N(C)Cc1ccccc1",
        "triethylamine": "CCN(CC)CC",
        "N,N-dimethylethanolamine": "CN(C)CCO",
        "invalid": "not_a_smiles"
    }
    
    for name, smi in examples.items():
        result, reason = is_tertiary_amine(smi)
        print(f"Example '{name}': {result} - {reason}")