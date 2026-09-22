
from Bio import SeqIO
from Bio.Data import CodonTable
import sys

# Candida albicans uses genetic code 12.
table = CodonTable.unambiguous_dna_by_id[12]


def read_dna(filename):
    record = SeqIO.read(filename, "fasta")
    return str(record.seq).upper().replace("U", "T")


def amino_acid(codon):
    if codon in table.stop_codons:
        return "*"
    return table.forward_table[codon]


# Read both alleles.
if len(sys.argv) != 3:
    sys.exit("Usage: python compare_alleles.py allele_A.fasta allele_B.fasta")

A = read_dna(sys.argv[1])
B = read_dna(sys.argv[2])

if len(A) != len(B):
    sys.exit("Sequences have different lengths. Align them first.")

if len(A) % 3 != 0:
    sys.exit("Sequence length is not divisible by 3.")

if set(A + B) - set("ACGT"):
    sys.exit("Sequences contain invalid or ambiguous nucleotides.")

# Compare the DNA sequences.
differences = []
synonymous = 0
missense = 0
other = 0

print("\nDNA DIFFERENCES")
print("-" * 65)

for i in range(len(A)):

    if A[i] == B[i]:
        continue

    # Find the codon containing this nucleotide.
    start = (i // 3) * 3

    codon_A = A[start:start + 3]
    codon_B = B[start:start + 3]

    aa_A = amino_acid(codon_A)
    aa_B = amino_acid(codon_B)

    aa_position = (i // 3) + 1

    # Determine the mutation type.
    if aa_A == aa_B:
        mutation = "Synonymous"
        synonymous += 1
    elif "*" in (aa_A, aa_B):
        mutation = "Stop codon change"
        other += 1
    else:
        mutation = "Missense"
        missense += 1

    differences.append(i + 1)

    print(
        f"Position {i + 1}: {A[i]} -> {B[i]} | "
        f"{codon_A} -> {codon_B} | "
        f"{aa_A}{aa_position}{aa_B} | {mutation}"
    )

# Translate both sequences into proteins.
protein_A = ""
protein_B = ""

for i in range(0, len(A), 3):
    protein_A += amino_acid(A[i:i + 3])
    protein_B += amino_acid(B[i:i + 3])

print("\nPROTEIN DIFFERENCES")
print("-" * 65)

protein_changes = 0

for i, (aa_A, aa_B) in enumerate(zip(protein_A, protein_B), 1):
    if aa_A != aa_B:
        print(f"Position {i}: {aa_A} -> {aa_B}")
        protein_changes += 1

if protein_changes == 0:
    print("No amino acid differences.")

# Print summary.
print("\nSUMMARY")
print("-" * 65)
print("Total nucleotide differences:", len(differences))
print("Synonymous substitutions:", synonymous)
print("Missense substitutions:", missense)
print("Stop codon changes:", other)
print("Amino acid differences:", protein_changes)
