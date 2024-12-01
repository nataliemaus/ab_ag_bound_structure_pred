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
    
    full_sequence = []
    chain_to_seq = {}
    
    # Iterate over all residues in the structure
    for model in structure:
        for chain in model:
            chain_to_seq[chain] = []
            for residue in chain:
                # Get the residue name (e.g., 'ALA', 'GLY', etc.)
                residue_name = residue.get_resname()
                
                # Map the three-letter residue name to one-letter code
                if residue_name in three_to_one:
                    full_sequence.append(three_to_one[residue_name])
                    chain_to_seq[chain].append(three_to_one[residue_name])
            chain_to_seq[chain] = ''.join(chain_to_seq[chain]) 
    
    # Join the sequence list into a string
    full_sequence = ''.join(full_sequence)

    return full_sequence, chain_to_seq

if __name__ == "__main__":
    # Replace with the path to your PDB file
    # pdb_file = "../her2/known_binding_pose/ag_only.pdb"
    # pdb_file = "../her2/known_binding_pose/ab_only.pdb"

    # pdb_file = "../influenza/known_pose_CR9114/inf_ag_only.pdb"
    pdb_file = "../influenza/known_pose_CR9114/inf_ab_only.pdb"

    # Extract and print the amino acid sequence
    amino_acid_sequence, chain_to_seq_dict = extract_sequence(pdb_file)
    print("Amino acid sequence:", amino_acid_sequence)
    print("chian to seq:")
    print(chain_to_seq_dict)

    
