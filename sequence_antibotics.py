# Autor：Evans
# Created time： 2021/11/11 下午 7:45

import copy

rna_amino_acid_code = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

dna_amino_acid_code = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                       'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                       'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                       'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                       'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                       'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                       'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                       'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

amino_acid_rna_code = {'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
                       'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'Y': ['UAU', 'UAC'],
                       '*': ['UAA', 'UAG', 'UGA'], 'C': ['UGU', 'UGC'], 'W': ['UGG'],
                       'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'],
                       'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'I': ['AUU', 'AUC', 'AUA'],
                       'M': ['AUG'], 'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'N': ['AAU', 'AAC'],
                       'K': ['AAA', 'AAG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'],
                       'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG']}

aa_table_dna = {'K': ['AAA', 'AAG'], 'N': ['AAC', 'AAT'], 'T': ['ACA', 'ACC', 'ACG', 'ACT'],
                'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
                'I': ['ATA', 'ATC', 'ATT'], 'M': ['ATG'], 'Q': ['CAA', 'CAG'], 'H': ['CAC', 'CAT'],
                'P': ['CCA', 'CCC', 'CCG', 'CCT'], 'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'E': ['GAA', 'GAG'],
                'D': ['GAC', 'GAT'], 'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'G': ['GGA', 'GGC', 'GGG', 'GGT'],
                'V': ['GTA', 'GTC', 'GTG', 'GTT'], '*': ['TAA', 'TAG', 'TGA'], 'Y': ['TAC', 'TAT'], 'C': ['TGC', 'TGT'],
                'W': ['TGG'], 'F': ['TTC', 'TTT']}

mass_table = {'G': '57', 'A': '71', 'S': '87', 'P': '97', 'V': '99', 'T': '101', 'C': '103', 'I': '113', 'L': '113',
              'N': '114', 'D': '115', 'K': '128', 'Q': '128', 'E': '129', 'M': '131', 'H': '137', 'F': '147',
              'R': '156', 'Y': '163', 'W': '186'}

def protein_translation(rna):
    rna_amino_acid_code = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
                           "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                           "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
                           "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
                           "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
                           "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                           "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                           "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                           "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
                           "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                           "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                           "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                           "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
                           "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                           "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                           "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
    protein_seq = ''
    for i in range(0,len(rna)-2,3):
        codon = rna[i:i+3]
        if rna_amino_acid_code[codon]=='*':
            break
        else:
            protein_seq += rna_amino_acid_code[codon]
    return protein_seq

def possibale_rna(peptide):
    num = 1
    for i in peptide:
        num *= len(amino_acid_rna_code[i])
    return num

def reversecompl(text):
    compl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    compl_txt = ''
    for i in text:
        compl_txt += compl[i]
    return compl_txt[::-1]

def dna_to_rna(dna):
    rna = ''
    for i in dna:
        if i == 'T':
            i = 'U'
            rna +=i
        else:
            rna +=i
    return rna

def peptide_encoding(string,peptide):
    dna_amino_acid_code = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S',
                           'TCG': 'S',
                           'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*',
                           'TGG': 'W',
                           'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
                           'CCG': 'P',
                           'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
                           'CGG': 'R',
                           'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
                           'ACG': 'T',
                           'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
                           'AGG': 'R',
                           'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
                           'GCG': 'A',
                           'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
                           'GGG': 'G'}
    peptide_len = len(peptide)
    peptide_substring = []
    reversed_string = reversecompl(string)
    #string with 3 read frames
    for j in range(3):
        protein_seq = ''
        if j == 0:
            j1 = 3
        else:
            j1 = j
        for i in range(j,len(string)-2-(3-j1),3):
            codon = string[i:i + 3]
            protein_seq += dna_amino_acid_code[codon]
        for i in range(len(protein_seq) - peptide_len + 1):
            peptides = protein_seq[i:i + peptide_len]
            if peptides == peptide:
                peptide_substring.append(string[j+i * 3:j+(i + peptide_len) * 3])
    #reverse string with 3 read frames
    for j in range(3):
        protein_seq = ''
        if j == 0:
            j1 = 3
        else:
            j1 = j
        for i in range(j,len(reversed_string)-2-(3-j1),3):
            codon = reversed_string[i:i + 3]
            protein_seq += dna_amino_acid_code[codon]
        for i in range(len(protein_seq) - peptide_len + 1):
            peptides = protein_seq[i:i + peptide_len]
            if peptides == peptide:
                peptide_substring.append(reversecompl(reversed_string[j+i * 3:j+(i + peptide_len) * 3]))
    return peptide_substring

def count_subpeptides(n):
    if n == 1:
        return 1
    else:
        result = n*(n-1)
    return result

def cycloSpectrum(peptide):
    mass_table = {'G': '57', 'A': '71', 'S': '87', 'P': '97', 'V': '99', 'T': '101', 'C': '103', 'I': '113', 'L': '113',
                  'N': '114', 'D': '115', 'K': '128', 'Q': '128', 'E': '129', 'M': '131', 'H': '137', 'F': '147',
                  'R': '156', 'Y': '163', 'W': '186'}
    mass_list = []
    for i in range(len(peptide)):
        #rotate the peptide
        rotate_peptide = peptide[i:] + peptide[:i]
        prefixMass = [0]
        for j in range(len(peptide)):
            prefixMass.append(int(prefixMass[j])+int(mass_table[rotate_peptide[j]]))
        total_mass = prefixMass[len(peptide)]
        for k in range(1,len(prefixMass)-1):
            mass_list.append(prefixMass[k])
    mass_list.append(0), mass_list.append(total_mass), mass_list.sort()
    return mass_list

def parentMass(peptide_ls):
    mass_list = []
    total_mass = sum(map(int,peptide_ls))
    for i in range(1,len(peptide_ls)):
        for j in range(len(peptide_ls)-i+1):
            peptide = peptide_ls[j:i+j]
            mass_list.append(sum(map(int,peptide)))
    mass_list.append(0), mass_list.append(total_mass), mass_list.sort()
    return mass_list

def CountingMass(Mass, masslist):
    Alphabet = {57: 'G', 71: 'A', 87: 'S', 97: 'P',
                99: 'V', 101: 'T', 103: 'C', 113: 'I/L',
                114: 'N', 115: 'D', 128: 'K/Q', 129: 'E',
                131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}
    if Mass == 0: return 1, masslist
    if Mass < 57: return 0, masslist
    if Mass in masslist: return masslist[Mass], masslist
    n = 0
    for i in Alphabet:
        k, masslist = CountingMass(Mass - i, masslist)
        n += k
    masslist[Mass] = n
    return n, masslist

def gauss_algorithm(n):
    return n*(n+1)/2

def is_sublist(sublist, list):
    if len(sublist) > len(list):
        return False
    for i in sublist:
        if i not in list:
            return False
        if sublist.count(i) > list.count(i):
            return False
    return True

#Code Challenge: Implement CyclopeptideSequencing (pseudocode reproduced below).
def CyclopeptideSequencing(Spectrum):
    mass_table = {57: 'G', 71: 'A', 87: 'S', 97: 'P',
                99: 'V', 101: 'T', 103: 'C', 113: 'I/L',
                114: 'N', 115: 'D', 128: 'K/Q', 129: 'E',
                131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}
    CandidatePeptides = {}
    final_peptides = set()
    Spectrum = list(map(str, Spectrum))
    #create a candidate peptide list which only includes the amnio acids in the spectrum
    #this list is also the base for candidate peptides to expand
    for i in Spectrum:
        if int(i) in mass_table:
            if i not in CandidatePeptides:
                CandidatePeptides[i] = [[i]]
            else:
                CandidatePeptides[i].append([i])
    # create a copy of CandidatePeptides
    # the copy is used to store the peptide with 1len which could used to expand the peptides
    CandidatePeptides_copy = copy.deepcopy(CandidatePeptides)
    def Expand(CandidatePeptides):
        CandidatePeptides_new = {}
        for key,values in CandidatePeptides.items():
            for value in values:
                for i in CandidatePeptides_copy:
                    new_key = int(key)+int(i)
                    #check whether the new_key(the mass of peptides) is in the spectrum
                    if str(new_key) in Spectrum:
                        #check whether the amino acide is about exceed the amount which spectrum have
                        if value.count(i) < len(CandidatePeptides_copy[i]):
                            new_value = copy.deepcopy(value)
                            new_value.append(i)
                            #use parentMass to generate the liner list of subpeptides mass
                            value_mass_str_ls = [str(i) for i in parentMass(new_value)]
                            #check if the mass of new peptide is consistent with the spectrum
                            if is_sublist(value_mass_str_ls,Spectrum):
                                if new_key not in CandidatePeptides_new:
                                    CandidatePeptides_new[new_key] = [new_value]
                                else:
                                    CandidatePeptides_new[new_key].append(new_value)
        return CandidatePeptides_new
    # expand the peptides till the CandidatePeptides only have peptides which mass is equal to total mass
    while len(CandidatePeptides)>1:
        CandidatePeptides = Expand(CandidatePeptides)
    # turn the peptides into strings result
    for values in CandidatePeptides.values():
        for value in values:
            final_peptides.add('-'.join(x for x in value))
    return final_peptides

def linearSpectrum(peptide):
    mass_table = {'G': '57', 'A': '71', 'S': '87', 'P': '97', 'V': '99', 'T': '101', 'C': '103', 'I': '113', 'L': '113',
                  'N': '114', 'D': '115', 'K': '128', 'Q': '128', 'E': '129', 'M': '131', 'H': '137', 'F': '147',
                  'R': '156', 'Y': '163', 'W': '186'}
    peptide_mass = [mass_table[i] for i in peptide]
    return parentMass(peptide_mass)




if __name__=='__main__':
    '''    
    rna = open('dataset_96_4.txt').read()
    data = open('result.txt','w')
    print(protein_translation(rna),file=data)
    data.close()
###
    handle = open('dataset_96_7.txt').read()
    handle = handle.split()
    string = handle[0]
    peptide = handle[1]
    data = open('result.txt','w')
    print('\n'.join(x for x in peptide_encoding(string,peptide)),file=data)
    data.close()
###
    print(count_subpeptides(21536))
###
    peptide = ''
    with open('dataset_98_4.txt') as handle:
        for line in handle:
            line = line.strip()
            peptide += line
    data = open('result.txt','w')
    print(*cycloSpectrum(peptide),file=data)
    data.close()
###
    peptide = ''print(CountingMass(1366, {})[0])
###
    n = 10379
    print(int(gauss_algorithm(n)+1))
###
    data = open('dataset_100_6.txt').read()
    Spectrum = []
    for i in data.split():
        Spectrum.append(i)
    data = open('result.txt','w')
    print(*CyclopeptideSequencing(Spectrum),file=data)
    data.close()
###
    peptide = 'HQPMLVDTHRATPKMHKLCRGQTMCLGMHEMWCYADSLTTRLLF'
    print(*linearSpectrum(peptide))
###
    ##quiz##
    string = ['CCUCGUACUGAUAUUAAU',
    'CCCAGUACCGAGAUGAAU',
    'CCCAGGACUGAGAUCAAU',
    'CCUCGUACAGAAAUCAAC']
    print('q2')
    for i in string:
        if protein_translation(i) == 'PRTEIN':
            bo = True
        else:
            bo = False
        print(i,protein_translation(i),bo)
    
    q3 = possibale_rna('LEADER')
    print('q3')
    print(q3)
    
    string = ['TAIM',
    'TMIA',
    'TLAM',
    'MTAI',
    'TALM',
    'MAIT'
    ]
    print('q5')
    for i in string:
        mass = ' '.join(str(x) for x in cycloSpectrum(i))
        if mass == '0 71 101 113 131 184 202 214 232 285 303 315 345 416':
            bo = True
        else:
            bo = False
        print(i,' '.join(str(x) for x in cycloSpectrum(i)),bo)
    
    print('q6')
    string = ['ETC',
    'CTQ',
    'AQV',
    'TCE',
    'CTV',
    'QCV']
    Spectrum_s = '0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333'
    # transfer Spectrum_s into a list
    Spectrum = [int(i) for i in Spectrum_s.split()]
    for i in string:
        mass = linearSpectrum(i)
        if is_sublist(mass,Spectrum):
            bo = True
        else:
            bo = False
        print(i,bo,mass)
    '''
