# Importing the necessary library for parsing the PDB file
from Bio.PDB import PDBParser

# Dictionary to map three-letter amino acid codes to one-letter codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def extract_sequence(pdb_file):
    # Create a PDB parser
    parser = PDBParser(QUIET=True)
    
    # Parse the structure from the PDB file
    structure = parser.get_structure("protein", pdb_file)
    
    sequence = []
    
    # Iterate over all residues in the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                # Get the residue name (e.g., 'ALA', 'GLY', etc.)
                residue_name = residue.get_resname()
                
                # Map the three-letter residue name to one-letter code
                if residue_name in three_to_one:
                    sequence.append(three_to_one[residue_name])
    
    # Join the sequence list into a string
    return ''.join(sequence)

if __name__ == "__main__":
    # Replace with the path to your PDB file
    pdb_file = "known_binding_pose/ag_only.pdb"
    # pdb_file = "known_binding_pose/ab_only.pdb"

    # Extract and print the amino acid sequence
    amino_acid_sequence = extract_sequence(pdb_file)
    print("Amino acid sequence:", amino_acid_sequence)
