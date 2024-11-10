#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
...Author: Shahid HADI...
This is the business logic API
"""

# Add the bl sub-directory to the module path (for testing this routine)
# and the directory above to import the config file
import sys
sys.path.insert(0, "../db/")
sys.path.insert(0, "../")

import dbapi   # Import the database api
import config  # Import configuration information (if needed)


# In[ ]:


def getAllEntries():
    """
    ...Function comment header goes here...
    This is a very simple function that just calls the database API to do the SQL to
    obtain the full list of entries. It doesn't need to do anything else.
    """
    return(dbapi.getAllEntries())


# In[ ]:


def searchGeneID(geneID):
    """
    ...Fnction to return query results for searching the database using gene ID as input...
    This function returns the 'searchEntryByGeneID(geneID)' from the database layer which takes
    its input from the serach box, provided by the cgi layer, for searching the database by using gene ID.
    
    It will return a list of list for all the matches as following:
        [("CTSL2","AB065492","mitochondrial ribosomal protein L41","9q33","MDPGLQQAL","aggagcttat agtctgaaaa")]
        
    If nothing is found it will return an empty list of lists
    """
    return(dbapi.searchEntryByGeneID(geneID))


# In[ ]:


def searchProteinProduct(proteinProductString):
    """
    ...Fnction to return query results for searching the database using protein product as input...
    This function returns the 'searchByProteinProduct(proteinProductString)' from the database layer which takes
    its input from the serach box, provided by the cgi layer, for searching the database by
    using protein product name.
    
    It will return a list of list for all the matches as following:
        [("CTSL2","AB065492","mitochondrial ribosomal protein L41","9q33","MDPGLQQAL","aggagcttat agtctgaaaa")]
        
    If nothing is found it will return an empty list of lists
    """
    return(dbapi.searchByProteinProduct(proteinProductString))


# In[ ]:


def searchAccession(accessionCode):
    """
    ...Fnction to return query results for searching the database using accession code as input...
    This function returns the 'searchByAccession(accessionCode)' from the database layer which takes
    its input from the serach box, provided by the cgi layer, for searching the database by
    using gene accession code.
    
    It will return a list of list for all the matches as following:
        [("CTSL2","AB065492","mitochondrial ribosomal protein L41","9q33","MDPGLQQAL","aggagcttat agtctgaaaa")]
        
    If nothing is found it will return an empty list of lists
    """
    return(dbapi.searchByAccession(accessionCode))


# In[ ]:


def searchChromLocation(chromosomeLocation):
    """
    ...Fnction to return query results for searching the database using  chromosome location as input...
    This function returns the 'searchByChromLocation(chromosomeLocation)' from the database layer which takes
    its input from the serach box, provided by the cgi layer, for searching the database by
    using chromosome location.
    
    It will return a list of list for all the matches as following:
        [("CTSL2","AB065492","mitochondrial ribosomal protein L41","9q33","MDPGLQQAL","aggagcttat agtctgaaaa")]
        
    If nothing is found it will return an empty list of lists
    """
    return(dbapi.searchByChromLocation(chromosomeLocation))


# In[ ]:


def getGeneralCodonFreq():
    """...Function to return overall chromosome codon frequencies, by apply the code in our 'geneCodonUsage()' function
    to the database...
    As this fuctions returns a set of fixed values, and for the purpose of speed and integrity, we applied our code
    to the database, and the returned values would be returned directly by this function without any calculations.
    
    Returns a dictionary with the keys being codes of the amino acids, like 'Ala(A)', each value
        consists of three lists:
            1 - List of codons used by the key amino acid.
            2 - List containing corresponding percentage usage for each codon.
            3 - List containing corresponding frequencies for each codon.     """
    return ({'Ala(A)': [('gct', 'gcc', 'gca', 'gcg'),
            [0.5, 1.35, 2.82, 1.75],
            [0.08, 0.21, 0.44, 0.27]],
            'Arg(R)': [('cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'),
            [1.07, 1.09, 0.62, 0.39, 0.86, 0.23],
            [0.25, 0.26, 0.15, 0.09, 0.2, 0.05]],
            'Asn(N)': [('aat', 'aac'), [2.13, 1.94], [0.52, 0.48]],
            'Asp(D)': [('gat', 'gac'), [2.11, 1.36], [0.61, 0.39]],
            'Cys(C)': [('tgt', 'tgc'), [1.37, 1.51], [0.48, 0.52]],
            'Gln(Q)': [('caa', 'cag'), [2.56, 0.95], [0.73, 0.27]],
            'Glu(E)': [('gaa', 'gag'), [2.48, 1.9], [0.57, 0.43]],
            'Gly(G)': [('ggt', 'ggc', 'gga', 'ggg'),
            [1.25, 1.12, 1.89, 0.95],
            [0.24, 0.22, 0.36, 0.18]],
            'His(H)': [('cat', 'cac'), [1.79, 0.93], [0.66, 0.34]],
            'Ile(I)': [('att', 'atc', 'ata'), [0.76, 3.41, 2.04], [0.12, 0.55, 0.33]],
            'Met(M)': ['atg', [4.54], [0.35]],
            'Leu(L)': [('tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'),
            [1.14, 3.09, 1.49, 1.8, 1.08, 2.35],
            [0.09, 0.24, 0.11, 0.14, 0.08, 0.55]],
            'Lys(K)': [('aaa', 'aag'), [1.89, 3.35], [0.45, 1.0]],
            'Phe(F)': [('ttt', 'ttc'), [3.37, 2.33], [0.59, 0.41]],
            'Pro(P)': [('cct', 'ccc', 'cca', 'ccg'),
            [0.48, 1.06, 1.7, 1.41],
            [0.1, 0.23, 0.37, 0.3]],
            'Ser(S)': [('tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'),
            [1.9, 1.31, 0.32, 1.27, 2.16, 1.86],
            [0.22, 0.15, 0.04, 0.14, 0.25, 0.21]],
            'Thr(T)': [('act', 'acc', 'aca', 'acg'),
            [0.46, 1.58, 2.07, 1.35],
            [0.08, 0.29, 0.38, 0.25]],
            'Trp(W)': ['tgg', [1.04], [1.0]],
            'Tyr(Y)': [('tat', 'tag'), [1.86, 1.66], [0.53, 0.47]],
            'Val(V)': [('gtt', 'gtc', 'gta', 'gtg'),
            [3.16, 0.86, 1.51, 1.11],
            [0.48, 0.13, 0.23, 0.17]],
            'STOP': [('taa', 'tga', 'tag'), [0.2, 0.09, 0.02], [0.66, 0.28, 0.07]]})


def getDNAseqAndCodingRegions(accession):
    """
    ...Function to return the coding regions, exons, represented within the entire gene DNA sequence...
    Input: Aceesion code of the gene of interest.
    
    This function generates the coding regions within a certain gene.
    
    Returns a list with multiple entries:
        The entries are sequences, with the same order of the gene.
        Entries that are simple strings, correspond to non coding regions and introns.
        Entries which are as strings in a list, within the original list, correspond to coding exons.
    """
# Assigning the results of the function of the db layer that  retruns all the details of a gene that
# BL laer needs to do its calculations
    gene = dbapi.getAllGeneInfo(accession)
        
# Putting the DNA gene sequence, cds list and protein amino acid sequence in to variables,
# calling them by there corresponding indices in result of the db layer function
    gene_seq = gene[1][0]
    cds = gene[1][1]
    protein_seq = gene[1][2]
        
# Substract 1 from each index number, indexing starts at zero in python, so, interpreting the indices given in the human
# readable file to python indexes, taking into account that when a string is sliced, the last element in the index is not
# returned!
    for i in range(0, len(cds), 2):
        cds[i] = cds[i] -1

# Adapting the cds to create a list of all seqs within the gene: non coding 5', cds 'its its exons and introns',
# non coding and 3' by creating a new list, complete_gene_seq_indices which holds the indices
# for every poisition within the gene.
    a = [0]
    if cds[0] != 0:
        complete_gene_seq_indices = a+cds
    else:
        complete_gene_seq_indices = cds
        
    if cds[-1] != len(gene_seq):
        complete_gene_seq_indices.append(len(gene_seq))
        
# Generatin all the subsequences that compose the gene:
    gene_seq_list = [gene_seq[i:j] for i,j in zip(complete_gene_seq_indices, complete_gene_seq_indices[1:])]

    


# Put the coding sequences in a list within the gene sequence list, upon requenst from front end!
    final_coding_seq_in_gene_seq = []
    for i, seq in enumerate(gene_seq_list):
        if cds[0] != 0:
            if i == 0 or i % 2 == 0:
                final_coding_seq_in_gene_seq.append(seq)
            elif i % 2 != 0:
                coding_seq_list = []
                coding_seq_list.append(seq)
                final_coding_seq_in_gene_seq.append(coding_seq_list)
        elif cds[0] == 0:
            if i == 0 or i % 2 == 0:
                coding_seq_list = []
                coding_seq_list.append(seq)
                final_coding_seq_in_gene_seq.append(coding_seq_list)
            elif i % 2 != 0:
                final_coding_seq_in_gene_seq.append(seq)
    return (final_coding_seq_in_gene_seq)


# In[ ]:


# Setting global variables for amino acid lists and corresponding codons list of lists that will be used
# by afterwards functions
single_letter_aa_list = ['A','R','N','D','C','Q','E','G','H','I','M','L','K','F','P','S','T','W','Y','V','STOP']

codon_list = [('gct', 'gcc', 'gca', 'gcg'),
 ('cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'),
 ('aat', 'aac'),
 ('gat', 'gac'),
 ('tgt', 'tgc'),
 ('caa', 'cag'),
 ('gaa', 'gag'),
 ('ggt', 'ggc', 'gga', 'ggg'),
 ('cat', 'cac'),
 ('att', 'atc', 'ata'),
 ('atg',),
 ('tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'),
 ('aaa', 'aag'),
 ('ttt', 'ttc'),
 ('cct', 'ccc', 'cca', 'ccg'),
 ('tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'),
 ('act', 'acc', 'aca', 'acg'),
 ('tgg',),
 ('tat', 'tag'),
 ('gtt', 'gtc', 'gta', 'gtg'),
 ('taa', 'tga', 'tag')]
# Creating a dictionary of amino acids and list of codons coding for any amino acid
aa_to_codon_dict = dict(zip(single_letter_aa_list, codon_list))


# In[ ]:


def alignNucleotide_to_aminoacid(accession):
    """
    ...Function to aligin coding DNA sequence to the amino acid sequnece of it protein product...
    Input: Accesion code of the gene of interest.
    
    This Function aligns the nucleotide sequence of a gene to the amino acid sequence of its protein.
    
    Retruns a list with each entry composed of the single letter amino acid and the codon being used to that
    amino acid, separated by a colon
    """
# Assigning the results of the function of the db layer that  retruns all the details of a gene that
# BL laer needs to do its calculations
    gene = dbapi.getAllGeneInfo(accession)
        
# Putting the DNA gene sequence, cds list and protein amino acid sequence in to variables,
# calling them by there corresponding indices in result of the db layer function
    gene_seq = gene[1][0]
    cds = gene[1][1]
    protein_seq = gene[1][2]
        
# Substract 1 from each index number, indexing starts at zero in python, so, interpreting the indices given in the human
# readable file to python indexes, taking into account that when a string is sliced, the last element in the index is not
# returned!
    for i in range(0, len(cds), 2):
        cds[i] = cds[i] -1


        
# Generate coding sequence string, coding_seq_str:
    coding_seq_str = ''
    for i in range (0, len(cds),2):
        code = gene_seq[cds[i]:cds[i+1]]
        coding_seq_str = coding_seq_str + code

# Generating two lists for codons in the coding DNA seq and a list of amino acids corresponding to the coding sequences
# The two lists are tied through their indices

    codon_gen_list =[]
    for i in range(0, len(coding_seq_str), 3):
        codon = coding_seq_str[i:i+3]
        codon_gen_list.append(codon)


    aa_gen_list = []
    for i in protein_seq:
        aa_gen_list.append(i)

# Generating a list of the coding sequence with the protein sequence separated by ':'
    string = ''
    aa_to_codon_list =[]

    for i,n in enumerate (aa_gen_list):
        string += n + ':' + codon_gen_list[i]
        aa_to_codon_list.append(string)
        string = ''
        
    
    return (aa_to_codon_list)


# In[ ]:


def RE_sites_list_cat(accession):
    """
    ...Function to retrun the retriction enzyme sites for 5 enzymes; EcoRI, BamHI, BsuMI, HindIII, and StuI,
    and categorize the according to their cuttig sites and applicability fo gene cloning purposes...
    Input: Accesion code of the gene of interest.
    
    This function Identifies the Restriction Enzymes "RE" valid for a gene, identiies the regions that each
    certain enzyme can cut, and categorising the enzymes, whether the restriction site for a certain RE
    falls within a coding region.
   
    Returns two dictionaries:
        1- First dictionary classifies the enzymes according to their availability and applicability. it has three lists:
            1- 'available_at_5p': should contain all the enzymes that cut at 5' end and DOES NOT cut at
                the coding region.
            2- 'available_at_3p': should contain all the enzymes that cut at 3' end and DOES NOT cut at
                the coding region.
            3- not_applicable: should contain all the enzymes that EITHER cut within the coding reion OR
                does not have any RE sites within the DNA sequence of the gene of interest
        2- Second dictionary links each enzyme to all the sites that it matches within the gene of interest
            by the first index of the match
    """
# Assigning the results of the function of the db layer that  retruns all the details of a gene that
# BL laer needs to do its calculations
    gene = dbapi.getAllGeneInfo(accession)
        
# Putting the DNA gene sequence, cds list and protein amino acid sequence in to variables,
# calling them by there corresponding indices in result of the db layer function
    gene_seq = gene[1][0]
    cds = gene[1][1]
    protein_seq = gene[1][2]
        
# Substract 1 from each index number, indexing starts at zero in python, so, interpreting the indices given in the human
# readable file to python indexes, taking into account that when a string is sliced, the last element in the index is not
# returned!
    for i in range(0, len(cds), 2):
        cds[i] = cds[i] -1

# Generate non-coding 5' sequence, entire coding sequence with introns included, and non-coding 3' sequence

    seq_5_prime = gene_seq[:cds[0]]
    complete_coding_seq = gene_seq[cds[0]:cds[-1]]
    seq_3_prime = gene_seq [cds[-1]:]

# Defining restriction enzymes:
    EcoRI = ['gaattc','']
    BamHI = ['ggatcc','']
    BsuMI = ['ctcgag','']
    HindIII = ['aagctt', '']
    StuI = ['aggcct', '']

    restriction_enzyme_list = [EcoRI, BamHI, BsuMI, HindIII, StuI]

    enzyme_name_list = ['EcoRI', 'BamHI', 'BsuMI', 'HindIII', 'StuI']

# Generating the palindromic sequence of each enzyme
    for ls in restriction_enzyme_list:
        i = 0
        while i < (len(ls)-1):
            ls[i+1] = ls[i+1].join(reversed(ls[i]))
            i+=1

# Generating lists to hold the indices for each enzyme
    EcoRI_sites = []
    BamHI_sites = []
    BsuMI_sites = []
    HindIII_sites = []
    StuI_sites = []
    RE_sites_list = [EcoRI_sites, BamHI_sites, BsuMI_sites, HindIII_sites, StuI_sites]

# Generating the lists for the categories that each enzyme falls into
    available_at_5p = []
    not_applicable = []
    available_at_3p = []
    availability_list = [available_at_5p, available_at_3p, not_applicable]
    availability_name_list = ['available_at_5p', 'available_at_3p', 'not_applicable']


    for n, enzyme in enumerate(restriction_enzyme_list):
        for m, seq in enumerate (enzyme):
            for i in range(len(gene_seq) - len(seq)):
                if gene_seq[i:i+len(seq)] == seq:
                    RE_sites_list[n].append(i)
    i = 0
    while i < len(restriction_enzyme_list):
        if not RE_sites_list[i]:
            availability_list[2].append(enzyme_name_list[i])
            i += 1
        else:
            if all(cds[0] < n >= cds[-1] for n in RE_sites_list[i]) and n < cds[0]:
                availability_list[0].append(enzyme_name_list[i])
    
            elif all(cds[0] < n >= cds[-1] for n in RE_sites_list[i]) and n > cds[-1]:
                availability_list[1].append(enzyme_name_list[i])
    
            else:
                availability_list[2].append(enzyme_name_list[i])
            i += 1



# Generating a list for the two dictioaries to be returned by the function
    RE_availability_dict = dict(zip(availability_name_list, availability_list))

    RE_sites_by_indices = dict(zip(enzyme_name_list, RE_sites_list))
    
    results_list = [RE_availability_dict, RE_sites_by_indices]
        
    return (results_list)


# In[ ]:


def geneCodonUsage(accession):
    """
    ...Function to return the codon percentage usage and codon frequecies for a gene of interest and
    to show statistical difference between the percentage codon usage of the entire chromosome and the
    gene of interest...
    Input: Accesion code of the gene of interest.
    
    This function calculates the codon percentage usage, ratio of a codon usage per 100 codons, and
    codon frequencies, the ratio of each codon amon codons used by a certain amino acid, and returns a list
    to show the difference between percentage codon usage of the entire chromosome 9 and the gene of interest
    
    Returns a dictionary with the keys being codes of the amino acids, like 'Ala(A)', each value
    consists of three lists:
        1 - List of codons used by the key amino acid.
        2 - List containing corresponding percentage usage for each codon.
        3 - List containing corresponding frequencies for each codon.
    And a list that holds the values of the statistical test for each corresponding codon.
    """
# Assigning the results of the function of the db layer that  retruns all the details of a gene that
# BL laer needs to do its calculations
    gene = dbapi.getAllGeneInfo(accession)
        
# Putting the DNA gene sequence, cds list and protein amino acid sequence in to variables,
# calling them by there corresponding indices in result of the db layer function
    gene_seq = gene[1][0]
    cds = gene[1][1]
    protein_seq = gene[1][2]
        
# Substract 1 from each index number, indexing starts at zero in python, so, interpreting the indices given in the human
# readable file to python indexes, taking into account that when a string is sliced, the last element in the index is not
# returned!
    for i in range(0, len(cds), 2):
        cds[i] = cds[i] -1
# Generate coding sequence
    coding_seq_str = ''
    for i in range (0, len(cds),2):
        code = gene_seq[cds[i]:cds[i+1]]
        coding_seq_str = coding_seq_str + code
# Creating a single list of all 64 codons
    all_codons_in_one_list = []

    for aa_list in codon_list:
        for codon in aa_list:
            all_codons_in_one_list.append(codon)
        
# capturing the indices of for every set of codons coding for a certain amino acid, to set a
# link between the indices of the all_codons_in_one_list and single_letter_aa_list
    codons_list_indices = [4,6,2,2,2,2,2,4,2,3,1,6,2,2,4,6,4,1,2,4,3]

    aggregate_codons_list_indices = []
    x = 0
    for i,n in enumerate (codons_list_indices):
        if i == 0:
            x = n
        else:
            x+=n
        aggregate_codons_list_indices.append(x)
        
# create all zero all_codons_count_list

    all_codons_count_list = []
    x = 0
    while x < 64:
        all_codons_count_list.append(0)
        x+=1

# Count every codons instances in the given gene

    for i in range(0, len(coding_seq_str), 3):
        codon = coding_seq_str[i:i+3]
        for x, y in enumerate (all_codons_in_one_list):
            if codon in y:
                all_codons_count_list[x]+=1

# Calculating percentage codon usages
    codons_percentage_usage_list = []

    for i in all_codons_count_list:
        codon_percentage = round((i / sum(all_codons_count_list)*100),2)
        codons_percentage_usage_list.append(codon_percentage)
    
# Separate percentages according to amino acids
    a = 0
    b = 0
    c = 1

    codon_percentage_map = []
    while c <21:
        if a == 0:
            count = codons_percentage_usage_list[:aggregate_codons_list_indices[a]]
            codon_percentage_map.append(count)
            a += 1
        else:
            count = codons_percentage_usage_list[aggregate_codons_list_indices[b]:aggregate_codons_list_indices[c]]
            codon_percentage_map.append(count)
            b+=1
            c+=1


# Saparate counts according to amino acids

    a = 0
    b = 0
    c = 1

    codon_count_map = []
    while c <21:
        if a == 0:
            count = all_codons_count_list[:aggregate_codons_list_indices[a]]
            codon_count_map.append(count)
            a += 1
        else:
            count = all_codons_count_list[aggregate_codons_list_indices[b]:aggregate_codons_list_indices[c]]
            codon_count_map.append(count)
            b+=1
            c+=1


# Change counts into frequencies:

    codon_freq_map = codon_count_map

    for ls in codon_freq_map:
        sum_ls = sum(ls)
        for i, instance in enumerate(ls):
            ls[i] =round ((instance/sum_ls),2)


# Combine codons count percentage and frequencies in on list

    count_freq_combined_list = []

    for i in range (len(codon_freq_map)):
        temp = []
        temp.append(codon_list[i])
        temp.append(codon_percentage_map[i])
        temp.append(codon_freq_map[i])
        count_freq_combined_list.append(temp)

# Create the final dictionary to be returned
    three_and_single_letter_aa_code = ['Ala(A)','Arg(R)','Asn(N)','Asp(D)','Cys(C)','Gln(Q)','Glu(E)',
                                   'Gly(G)','His(H)','Ile(I)','Met(M)','Leu(L)','Lys(K)','Phe(F)',
                                   'Pro(P)','Ser(S)','Thr(T)','Trp(W)','Tyr(Y)','Val(V)','STOP']

    codon_usage_dict = dict(zip(three_and_single_letter_aa_code, count_freq_combined_list))

    
# Setting a statistical test to show significant different between percentage codon usages of a single gene and the
# entire chromosome 9. the test returns one if the difference is greater than two fold and zero if smaller than two
# fold.

# Chrom_9_codon_percentage_map list holds the percentage codon usage values for the entire chromosome 9
    Chrom_9_codon_percentage_map = [[0.5, 1.35, 2.82, 1.75], [1.07, 1.09, 0.62, 0.39, 0.86, 0.23], [2.13, 1.94],
                                [2.11, 1.36], [1.37, 1.51], [2.56, 0.95], [2.48, 1.9], [1.25, 1.12, 1.89, 0.95],
                                [1.79, 0.93], [0.76, 3.41, 2.04], [4.54], [1.14, 3.09, 1.49, 1.8, 1.08, 2.35],
                                [1.89, 3.35], [3.37, 2.33], [0.48, 1.06, 1.7, 1.41],
                                [1.9, 1.31, 0.32, 1.27, 2.16, 1.86], [0.46, 1.58, 2.07, 1.35], [1.04],
                                [1.86, 1.66], [3.16, 0.86, 1.51, 1.11], [0.2, 0.09, 0.02]]

# stats_list created to hold the results returned by our statistical test
    stats_list = [[0,0,0,0], [0, 0, 0, 0, 0, 0], [0,0 ], [0,0], [0,0], [0,0], [0,0], [0,0,0,0], [0,0], [0,0,0], [0],
              [0,0 ,0 ,0 ,0 ,0 ], [0,0], [0,0], [0,0,0,0], [0,0,0,0,0,0], [0,0 ,0 ,0 ], [0], [0,0], [0,0 ,0 ,0 ],
              [0,0 ,0 ]]

# Statistical test desigend to show more than two fold difference between the percentage codon usage
# of entire chromosom 9 and the gene of interest
    for i, n in enumerate(Chrom_9_codon_percentage_map):
        for e, m in enumerate(n):
            if Chrom_9_codon_percentage_map[i][e] != 0 and codon_percentage_map[i][e] != 0:
                if 0.5 < Chrom_9_codon_percentage_map[i][e]/codon_percentage_map[i][e] < 2:
                    stats_list[i][e] = 0
                else:
                    stats_list[i][e] = 1
            elif Chrom_9_codon_percentage_map[i][e] == 0 or codon_percentage_map[i][e] == 0:
                stats_list[i][e] = 1

    return (codon_usage_dict, stats_list)


# In[1]:





# In[ ]:


def input_seq_RE_sites(accession, your_choice):
    """
    ...Function to retrun results for a restriction enzyme sequence entered by user...
    Input: accesion: accession code of gene of interest
            your_choice: sequence as a string entered by a user, obtained from front layer
            
     Returns two dictionaries:
        1- First dictionary classifies the input sequence according to their availability and applicability.
           it has three lists:
            1- 'available_at_5p': should contain the user's choice if it cuts at 5' end and DOES NOT cut at
                the coding region.
            2- 'available_at_3p': should contain the user's choice cut at 3' end and DOES NOT cut at
                the coding region.
            3- not_applicable: should contain the user's choice if it EITHER cut within the coding reion OR
                does not have any RE sites within the DNA sequence of the gene of interest
        2- Second dictionary links the user choice to all the sites that it matches within the gene of interest
            by the first index of the match
         """
# Assigning the results of the function of the db layer that  retruns all the details of a gene that
# BL laer needs to do its calculations
    gene = dbapi.getAllGeneInfo(accession)
        
# Putting the DNA gene sequence, cds list and protein amino acid sequence in to variables,
# calling them by there corresponding indices in result of the db layer function
    gene_seq = gene[1][0]
    cds = gene[1][1]
    protein_seq = gene[1][2]
        
# Substract 1 from each index number, indexing starts at zero in python, so, interpreting the indices given in the human
# readable file to python indexes, taking into account that when a string is sliced, the last element in the index is not
# returned!
    for i in range(0, len(cds), 2):
        cds[i] = cds[i] -1

# Commenting for rhis function is present at previous 'RE_sites_list_cat()' function!
    for i in range(0, len(cds), 2):
        cds[i] = cds[i] -1
    input_seq = your_choice.replace('-','').lower()
    

    seq_5_prime = gene_seq[:cds[0]]
    complete_coding_seq = gene_seq[cds[0]:cds[-1]]
    seq_3_prime = gene_seq [cds[-1]:]

    input_seq_list = [input_seq, '']
    input_seq_list[1] = input_seq_list[1].join(reversed(input_seq_list[0]))

    input_seq_sites = []

    for i in range(len(gene_seq) - len(input_seq)):
            if gene_seq[i:i+len(input_seq)] == input_seq:
                input_seq_sites.append(i)

    available_at_5p = []
    not_applicable = []
    available_at_3p = []
    availability_list = [available_at_5p, available_at_3p, not_applicable]
    availability_name_list = ['available_at_5p', 'available_at_3p', 'not_applicable']


    if not input_seq_sites:
        availability_list[2].append("input_seq")
    else:
        if all(cds[0] < n >= cds[-1] for n in input_seq_sites) and n < cds[0]:
            availability_list[0].append("your sequence")
    
        elif all(cds[0] < n >= cds[-1] for n in input_seq_sites) and n > cds[-1]:
            availability_list[1].append("your sequence")
    
        else:
            availability_list[2].append("your sequence")

    
    input_seq_site_dict = {"Your sequence sites": input_seq_sites}

    RE_availability_dict = dict(zip(availability_name_list, availability_list))
    
    results_list = [RE_availability_dict, input_seq_site_dict]

    return(results_list)
    
    


