import collections
from Utilities import codon_table

DNA_ReverseComplement = {"A": "T", "C": "G", "G": "C", "T": "A"}
Nucleotides = list(DNA_ReverseComplement.keys())

# Check the sequence to make sure it is a DNA String
def validateSeq (dna_seq) :
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
          return False
    return tmpseq

def countNucFrequency (dna_seq):
    return dict(collections.Counter(dna_seq))

def transcription (dna_seq):
    return dna_seq.replace("T", "U")

def reverseComplement (dna_seq):
    return ''.join([DNA_ReverseComplement[nuc] for nuc in dna_seq])[::-1]

def getGC(data, nuc1 = 'G', nuc2 = 'C'):
    allNuc = countNucFrequency(data)
    return (allNuc[nuc1] + allNuc[nuc2]) / sum(allNuc.values()) * 100

def hamming_distance(s, t):
    if len(s) != len(t):
        raise ValueError("Strands must be of equal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s, t))

def transform(dna_seq):
    return ''.join([codon_table[dna_seq[n:n+3]] for n in range(0, len(dna_seq), 3) if codon_table[dna_seq[n:n+3]] != 'Stop'])


def formular(total, a, b, multiplier = 1 ):
    return (a/total) * (b/(total-1)) * multiplier

def mendelian(k, m, n):
    total = k + m + n
    prob = 0
    prob += formular(total, k, k-1)
    prob += formular(total, k, m)
    prob += formular(total, k, n)
    prob += formular(total, m, k)
    prob += formular(total, m, m-1, 0.75)
    prob += formular(total, m, n, 0.5)
    prob += formular(total, n, k)
    prob += formular(total, n, m, 0.5)
    prob += formular(total, n, n-1, 0)
    return prob

def gestAllSubstrings(dna, sub):
    index = []
    for i in range(len(dna)):
        if dna[i:].startswith(sub):
            index.append(i+1)
    return index





def readFile(filePath):
    with open(filePath, 'r') as f:
        return [l.strip() for l in f.readlines()]
            

def getFastaDict(data, fastaDict = {}):  
    for line in data:
        if '>' in line:
            fastaDict[line] = ""
            currentLabel = line
        else:
            fastaDict[currentLabel] += line.strip()
    return fastaDict
 

def getConsensusString(FASTADict, profileMatrix = {}, consensusString = ""):
    for value in FASTADict.values():
        for position in range(len(value)):
            if position not in profileMatrix:
                profileMatrix[position] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            profileMatrix[position][value[position]] += 1

    for position in profileMatrix:
        maxCount = 0
        maxLetter = ""
        for letter in profileMatrix[position]:
            
            if profileMatrix[position][letter] > maxCount:
                maxCount = profileMatrix[position][letter]
                maxLetter = letter
        consensusString += maxLetter
    return consensusString, profileMatrix


def printProfileMatrix(profileMatrix):
    a,c,g,t = [],[],[],[]

    for letter in profileMatrix:
        a.append(profileMatrix[letter]['A'])
        c.append(profileMatrix[letter]['C'])
        g.append(profileMatrix[letter]['G'])
        t.append(profileMatrix[letter]['T'])
    return {"A":a, "C":c, "G":g, "T":t}
        
'''
FASTAFile = readFile( "Main/Input.txt")
FASTADict = getFastaDict(FASTAFile)
output = getConsensusString(FASTADict)
matrixToPrint = printProfileMatrix(output[1])
print(output[0] + "\n", *(f"{i}: {' '.join(map(str, matrixToPrint[i]))}\n" for i in matrixToPrint))
'''


def getOverlapGraph(FASTADict):
    for key1 in FASTADict:
        for key2 in FASTADict:
            if key1 != key2:
                if FASTADict[key1][-3:] == FASTADict[key2][:3]:
                    # remove > from key
                    print(key1[1:], key2[1:])
    return


def expectedOffspring(couples):
    return couples[0] * 2 + couples[1] * 2 + couples[2] * 2 + couples[3] * 1.5 + couples[4] * 1 + couples[5] * 0


def mendelsSecondLaw(k, N):
    total = 2 ** k
    dominant = 0
    for i in range(N, total + 1):
        dominant += binomial(total, i) * (0.25 ** i) * (0.75 ** (total - i))
    return dominant

def binomial(n, k):
    return factorial(n) / (factorial(k) * factorial(n - k))

def factorial(n):
    if n == 0:
        return 1
    return n * factorial(n - 1)