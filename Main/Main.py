from Toolkit import *
from Utilities import *
import random

# Creating a random DNA sequence for testing:
randDNAStr = ''.join( [random.choice(Nucleotides)
                       for nuc in range(50)])

DNAStr = validateSeq (randDNAStr)

print(f'\nSequence: {colored(DNAStr)}\n')
print(f'[1] + Sequence Length: {len(DNAStr)} \n ')
print(colored(f'[2] + Nucleotide Frequency: {countNucFrequency(DNAStr)}\n '))
print(f'[3] + DNA/RNA Transcription: {colored(transcription(DNAStr)) } \n ')
print(f"[4] + DNA String + Reverse Complement: \n5' {colored(DNAStr)} 3'")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {colored(reverseComplement(DNAStr))} 5'\n") 