# Find-Homopolymer-Indels
These programs can extend trinucleotide sequences and count homopolymer sequences. The homopolymer length is recorded in a seperate BED file, along with the homopolymer type. Homopolymers that have less than four of the same nucleotide in a row were decided to be non-homopolymers and marked as 'NA' on the homopolymer type.

## Requirements
- Cromosomes in txt files
- BED file in standard format with a starting trinucleotide sequence
- Pandas package

## Other
These programs were made to analyze UV mutated DNA from S. cerevisiae. Data was supplemented by the Wyrick Laboratory at Washington State University.
