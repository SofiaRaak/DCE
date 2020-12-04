from Bio import SeqIO
impurt numpy as np

#Dictionary determining molecular weight of each amino acid in kDa
fasta_to_weight = {'A':0.0891, 'R':0.1742, 'N': 0.1321, 'D':0.1331, 'C': 0.1212,
                  'E': 0.1471, 'Q': 0.1462, 'G': 0.0751, 'H': 0.1552, 'I': 0.1312,
                  'L': 0.1312, 'K': 0.1462, 'M': 0.1492, 'F': 0.1652, 'P': 0.1151,
                  'S': 0.1051, 'T': 0.1191, 'W': 0.2042, 'Y': 0.1812, 'V': 0.1171}

#Separate sequences, headers
def SeparateSequences(file, filetype):
    
    global ID_list
    global seq_list
    
    ID_list = []
    seq_list = []
    
    for seq_record in SeqIO.parse(file, filetype):
        ID_list.append(seq_record.id)
        seq_list.append(seq_record.seq)


#Find molecular weight
def FindkDa(fasta):
    
    mol_weights = []
    
    for letter in protein:
        letter = letter.upper()
        mol_kDa = fasta_to_weight[letter]
        mol_weights.append(mol_kDa)
    
    return round(sum(mol_weights))


#Output molecular weight for each protein in list
def MolWeights(file, filetype):
    
    SeparateSequences(file, filetype)
    
    mol_weights = np.zeros([len(seq_list), 2])
    
    for i in range(len(seq_list)):
        mol_weights[i, 0] = ID_list[i]
        mol_weights[i, 1] = 