import itertools

def translate_dna():
    seqm = input("Input DNA = ").upper()
    print(" ")

    codon_map = {
        'UUU': ['Phenylalanine', 'F', 'Phe'], 'UUC': ['Phenylalanine', 'F', 'Phe'], 'UUA': ['Leucine', 'L', 'Leu'], 'UUG': ['Leucine', 'L', 'Leu'], 
        'CUA': ['Leucine', 'L', 'Leu'], 'CUG': ['Leucine', 'L', 'Leu'], 'AUU': ['Isoleucine', 'I', 'Ile'], 'AUC': ['Isoleucine', 'I', 'Ile'], 
        'AUA': ['Isoleucine', 'I', 'Ile'], 'AUG': ['Methionine', 'M', 'Met'], 'GUU': ['Valine', 'V', 'Val'], 'GUC': ['Valine', 'V', 'Val'], 
        'GUA': ['Valine', 'V', 'Val'], 'GUG': ['Valine', 'V', 'Val'], 'UCU': ['Serine', 'S', 'Ser'], 'UCC': ['Serine', 'S', 'Ser'], 
        'UCA': ['Serine', 'S', 'Ser'], 'UCG': ['Serine', 'S', 'Ser'], 'CCU': ['Proline', 'P', 'Pro'], 'CCC': ['Proline', 'P', 'Pro'], 
        'CCA': ['Proline', 'P', 'Pro'], 'CCG': ['Proline', 'P', 'Pro'], 'ACU': ['Threonine', 'T', 'Thr'], 'ACC': ['Threonine', 'T', 'Thr'], 
        'ACA': ['Threonine', 'T', 'Thr'], 'ACG': ['Threonine', 'T', 'Thr'], 'GCU': ['Alanine', 'A', 'Ala'], 'GCC': ['Alanine', 'A', 'Ala'], 
        'GCA': ['Alanine', 'A', 'Ala'], 'GCG': ['Alanine', 'A', 'Ala'], 'UAU': ['Tyrosine', 'Y', 'Tyr'], 'UAC': ['Tyrosine', 'Y', 'Tyr'], 
        'UAA': ['Stop', 'end', 'Stop'], 'UAG': ['Stop', 'end', 'Stop'], 'UGA': ['Stop', 'end', 'Stop'], 'UGG': ['Tryptophan', 'W', 'Trp'], 
        'CAA': ['Glutamine', 'Q', 'Gln'], 'CAG': ['Glutamine', 'Q', 'Gln'], 'AAU': ['Asparagine', 'N', 'Asn'], 'AAC': ['Asparagine', 'N', 'Asn'], 
        'AAA': ['Lysine', 'K', 'Lys'], 'AAG': ['Lysine', 'K', 'Lys'], 'GAU': ['Aspartic Acid', 'D', 'Asp'], 'GAC': ['Aspartic Acid', 'D', 'Asp'], 
        'GAA': ['Glutamic Acid', 'E', 'Glu'], 'GAG': ['Glutamic Acid', 'E', 'Glu'], 'UGU': ['Cysteine', 'C', 'Cys'], 'UGC': ['Cysteine', 'C', 'Cys'], 
        'CGA': ['Arginine', 'R', 'Arg'], 'CGG': ['Arginine', 'R', 'Arg'], 'AGU': ['Serine', 'S', 'Ser'], 'AGC': ['Serine', 'S', 'Ser'], 
        'AGA': ['Arginine', 'R', 'Arg'], 'AGG': ['Arginine', 'R', 'Arg'], 'GGU': ['Glycine', 'G', 'Gly'], 'GGC': ['Glycine', 'G', 'Gly'], 
        'GGA': ['Glycine', 'G', 'Gly'], 'GGG': ['Glycine', 'G', 'Gly']
    }

    # Validate DNA sequence length
    if len(seqm) % 3 != 0:
        print("Nucleotides should be in multiples of 3.")
        print(" ")
        return

    # Validate if the input is a valid DNA sequence
    valid = seqm.count("A") + seqm.count("C") + seqm.count("T") + seqm.count("G")
    if valid == len(seqm):
        complement_seq = ""
        mrna_seq = ""
        for i in seqm:
            if i == "A":
                complement_seq += "T"
                mrna_seq += "U"
            elif i == "T":
                complement_seq += "A"
                mrna_seq += "A"
            elif i == "C":
                complement_seq += "G"
                mrna_seq += "G"
            elif i == "G":
                complement_seq += "C"
                mrna_seq += "C"

        print("Complement = ", complement_seq)
        print("mRNA = ", mrna_seq)

        amino_acids = []
        for j in range(0, len(mrna_seq) - 2, 3):
            codon = mrna_seq[j:j + 3]
            if codon in codon_map:
                if codon_map[codon][0] == 'Stop':
                    break
                amino_acids.append(codon_map[codon][2] + " (" + codon_map[codon][1] + ")")
        
        print("Aminoacids = ", " - ".join(amino_acids))
        print(" ")

    # Proceed to the next step: Get amino acid input and find codon frequencies
    get_mrna_and_codon_frequency(codon_map)

# Function to get the codon frequencies based on input amino acids
def get_mrna_and_codon_frequency(codon_map):
    amino_acid_to_codon = {}
    
    # Reverse map amino acid symbols (single letter) to their codons
    for codon, amino_info in codon_map.items():
        amino_acid_symbol = amino_info[1]
        if amino_acid_symbol not in amino_acid_to_codon:
            amino_acid_to_codon[amino_acid_symbol] = []
        amino_acid_to_codon[amino_acid_symbol].append(codon)

    # Get input amino acids
    input_aa = input("Input Aminoacid (max 3 letters, e.g. WYW) = ").upper()

    if len(input_aa) > 3:
        print("Error: You can only input a maximum of 3 amino acids.")
        return

    # Generate all possible mRNA sequences by using itertools.product to get all codon combinations
    codon_combinations = [amino_acid_to_codon[aa] for aa in input_aa if aa in amino_acid_to_codon]
    possible_mrna_sequences = list(itertools.product(*codon_combinations))

    for mRNA_seq in possible_mrna_sequences:
        joined_mrna_seq = "".join(mRNA_seq)
        print("\nmRNA = ", joined_mrna_seq)

        # Count the frequency of codons in the mRNA sequence
        codon_frequency = {}
        for codon in mRNA_seq:
            if codon in codon_frequency:
                codon_frequency[codon] += 1
            else:
                codon_frequency[codon] = 1

        # Print the frequency of each codon
        for codon, freq in codon_frequency.items():
            print(f"{codon} = {freq}")

# Run the main function
translate_dna()
