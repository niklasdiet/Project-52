from Toolkit import *
from Utilities import *
import requests

FASTAFile = readFile("Main/Input.txt")


for id in FASTAFile:
    url = f'http://www.uniprot.org/uniprot/{id}.fasta'
    #call the url and get the fasta file
    fasta = requests.get(url).text
    print(fasta)