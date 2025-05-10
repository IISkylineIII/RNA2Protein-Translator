# RNA2Protein-Translator

# Description
This script translates an RNA sequence into a protein sequence based on a given genetic code dictionary. It reads codons (triplets of nucleotides), maps them to their corresponding amino acids, and stops translation when a stop codon is encountered.

The function mimics the biological process of translation and is commonly used in computational biology pipelines for transcriptome and protein synthesis simulations.

# Usage
Example:

```
rna_string = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
genetic_code = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

amino_acid_string = translate_rna_to_protein(rna_string, genetic_code)
print(amino_acid_string)
```
# Output:
MAMAPRTEINSTRING

# Function
```
def translate_rna_to_protein(rna, genetic_code):
    protein = ""

    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        amino_acid = genetic_code[codon]

        if amino_acid == "STOP":
            break

        protein += amino_acid

    return protein

```

# Applications

* RNA-to-Protein translation in transcriptomics and functional genomics.
* Synthetic biology for designing expression constructs.
* Teaching bioinformatics and molecular biology.
* Used in pipeline stages for proteomic annotation.

# License
This project is licensed under the MIT License.










