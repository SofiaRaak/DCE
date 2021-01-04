from Bio import Entrez
from Bio import SeqIO
import pandas as pd

#Dictionary determining molecular weight of each amino acid in kDa
fasta_to_weight = {'A':0.0891, 'R':0.1742, 'N': 0.1321, 'D':0.1331, 'C': 0.1212,
                  'E': 0.1471, 'Q': 0.1462, 'G': 0.0751, 'H': 0.1552, 'I': 0.1312,
                  'L': 0.1312, 'K': 0.1462, 'M': 0.1492, 'F': 0.1652, 'P': 0.1151,
                  'S': 0.1051, 'T': 0.1191, 'W': 0.2042, 'Y': 0.1812, 'V': 0.1171}


#Separate sequences, headers
def SeparateSequences(file):
    
    global ID_list
    global seq_list
    
    ID_list = []
    seq_list = []
        
    for seq_record in SeqIO.parse(file, 'fasta'):
        ID_list.append(seq_record.id)
        seq_list.append(seq_record.seq)

#Download sequence files with Entrez
def DownloadSequence(my_email, acc_numbers): 

    Entrez.email = my_email
    handle = Entrez.efetch(db='protein', id=acc_numbers, rettype='fasta', retmode='text')

    SeparateSequences(handle)

    
#Find molecular weight
def FindkDa(protein):
    
    mol_weights = []
    
    for letter in protein:
        letter = letter.upper()
        mol_kDa = fasta_to_weight[letter]
        mol_weights.append(mol_kDa)
    
    #Likely overestimation due to excess water in the protein, but good enough for crude estimation
    #Overestimation likely minimal if post-translational modifications present
    return sum(mol_weights)


#Output molecular weight for each protein in list
def MolWeights(my_email, acc_numbers, download=True, file=None):
    
    if download == True:
        DownloadSequence(my_email, acc_numbers)
    
    elif download == False:
        SeparateSequences(file)
    
    mol_weights = pd.DataFrame(columns = ['ID', 'kDa'])
    
    for i in range(len(seq_list)):
        mol_weights.loc[i, 'ID'] = ID_list[i]
        mol_weights.loc[i, 'kDa'] = FindkDa(seq_list[i])
        
    return mol_weights
