#!/usr/bin/env python

import sys
import argparse

################################################################################################
################################################################################################

print("")
print("")
print("")

################################################################################################
################################################################################################ PART 1

# 1. Write a function that takes as an input a DNA string and returns the reverse complement of that string.

#    a. Write the function assuming you will only get ‘A’, ‘C’, ‘T’, and ‘G’ characters.
#    b. What other characters beside ‘A’, ‘C’, ‘T’, and ‘G’ can be found in fastq or fasta files?
#    c. What is the rationale for the use of non ‘A’, ‘C’, ‘T’, and ‘G’ characters? State some ways you might extend the reverse complement function to handle these characters.

################################################################################################
################################################################################################ PART 1.a

# a. Write the function assuming you will only get ‘A’, ‘C’, ‘T’, and ‘G’ characters.

# to define the function to REVERSE_COMPLEMENT the string
def reversecomplement(seq):
    '''This function returns a reverse complement of a DNA sequence'''

    # Complement the DNA strand
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()

    # Reverse the complemented DNA strand
    seq = seq[::-1]
    return seq

# to verify that the sequence was correctly reverse complemented :
# providing an example
# seq='ATCGATCAGTCCTAGCATCG'
# print(seq)
# print('compute the revese complement of a sequence :', seq)
# print('the reverse complement is :')
# print(reversecomplement(seq))
# print('whether the reverse complement is correct :')
# print(reversecomplement(seq) == 'CGATGCTAGGACTGATCGAT')

################################################################################################
################################################################################################ PART 1.b

# b. What other characters beside ‘A’, ‘C’, ‘T’, and ‘G’ can be found in fastq or fasta files?

# Let's consider as an example the following sequence :

# @105220b1-f48a-43b4-8e89-a3cf20afeb0d runid=fa1d76e661ac2bbb53a002e85e75a30e91827c51 sampleid=1 read=5559 ch=33 start_time=2019-10-18T22:15:11Z
# CAAGCATACTTCATTCAGTCAGGCGAAATTATTGCCAGGTCGCCGCCTACCGTGACAAGAAAGTTGTCCGGTATCTTTGTGTTTCTGTTGGTGCTGATATTGCGTTATGCATGAACGTAATGCTCATCGGATTGTGAATCCACCATGCGCGGAAAGGCAGGGCGACAGGCAAGTCACAAGGAACACCAGACGCTTTGT
# +
# ,&'%/..0013618'.'(*'#(%#$&&%%,-//&#'$&&'&%$('''4+9;7<87756-1-+$$,3%%;)-%$%%$&)''#067897:9$%)(<AC?,(**$&'.6:<=394./41*,12((:?;3/9'(-4<)=99D99BFAC9:;588+.&+&%()%-(59,,($,12+,,*6(-))&-$&0'%)%&%'""%(0((

# In a FASTA file, we may find the following characters :

# N: an ambiguous nucleotide (A, C, T, or G)
# IUPAC ambiguity codes:
#        R (A or G)
#        Y (C or T)
#        M (A or C)
#        K (G or T)
#        S (C or G)
#        W (A or T)
#        B (not A)
#        D (not C )
#        H (not G)
#        V (not T)
# - : Gap characters, indicating insertions or deletions in alignments.

# In a FASTQ file :

# In the field 1, @' character is followed by a sequence identifier (a string of ASCII characters)
# In the field 2, the sequence of the base calls; A, C, T, G and N (the ambiguous nucleotide).
# In the field 3, a plus (+) sign as separator.
# In the field 4, the Quality Score characters.
# The quality scores give an indication of the confidence level in each base call.
# The quality scores are ASCII characters ranging from '!' (lowest quality) to '~' (highest quality).
# The specific encoding, such as Phred+33 or Phred+64, helps interpret the quality scores.

# In a "Phred+33" specification (Sanger format), the ASCII value is offset by 33.
# The range of characters from '!' (ASCII 33) to 'I' (ASCII 73) represents quality scores from 0 to 40.

################################################################################################
################################################################################################ PART 1.c

# c. What is the rationale for the use of non ‘A’, ‘C’, ‘T’, and ‘G’ characters?
# State some ways you might extend the reverse complement function to handle these characters.

# Considering the case of "N" (representing ambiguity) or '-' (representing deletions or insertions),
# in a reverse complemented sequence, we could use :
# "N" as the reverse complement af "N", and "-" (gap) as a reverse complement character of "-".

# In this case c, we modify the reverse complement sequence in the following way :

def reverse_complement2(seq):
    '''This function returns a reverse complement of a DNA sequence'''

    # Complement the DNA strand
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").replace("N", "n").replace("-", "-")
    seq = seq.upper()

    # Reverse the complemented DNA strand
    seq = seq[::-1]
    return seq

# to verify that the sequence was correctly reverse complemented :
# providing an example :
# seq='ATNCGATCAGTCCTAGCAT-CG'
# print(seq)
# print('compute the revese complement of a sequence :', seq)
# print('the reverse complement is :')
# print(reverse_complement2(seq))
# print('whether the reverse complement is correct :')
# print(reverse_complement2(seq) == 'CG-ATGCTAGGACTGATCGNAT')

################################################################################################
################################################################################################ PART 2
# 2. Write a function which returns True if a DNA string is a palindrome or False if it is not a palindrome.

# In a palindrome sequence, the reverse complement of the sequence must read the same as the original sequence.
# For example, the sequence 5'-CGATCG-3' is considered a palindrome since its reverse complement 3'-GCTAGC-5' reads the same.

# We add another function to our code that check whether a sequence is a palindrome sequence :

def palindrome_sequence(seq):
    '''This function checks if the given DNA sequence is a palindrome'''

    # Get the reverse complement of the sequence
    rev_comp = reverse_complement2(seq)
    rev_comp = rev_comp.upper()

    # Check if the sequence is the same as its reverse complement
    return seq == rev_comp

# To verify the results :
# string = "GATTACA"
# print("To verify if the function palindrome_sequence() works correctly on the string:", string)
# print("Is the sequence a palindrome?", palindrome_sequence(string))

# print("To verify if the function palindrome_sequence() works correctly on the string:", string)
# print("Is the sequence a palindrome?", palindrome_sequence(string))

################################################################################################
################################################################################################ PART 3.a
# 3. Given the following DNA string: TGCAGGACGTATGGCCGTTAGCTAACGGCCATAGATGCTATG
#
#    a. Write a function which counts kmers for any given length (k). What is the most frequent 3, 4 and 5nt length kmer found in the sequence above?
#    b. Write a function which identifies the longest kmer in the string above that is also a palindrome.

sequence = "TGCAGGACGTATGGCCGTTAGCTAACGGCCATAGATGCTATG"

# The function to count K-mers is :

def count_kmers(sequence, k):
    kmer_counts = {}                       # an empty dictionary where we keep the counts for a kmer of length k
    for i in range(len(sequence) - k + 1): # we compute the length of the shortened sequence
        kmer = sequence[i:i+k]             # we extract a sub-string from position "i" to "position "i+k"
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1
    return kmer_counts

################################################################################################

print("Given the sequence in the HW :", sequence)
k = 3
print("the frequency of the k-mers of length", k)
print(count_kmers(sequence, k))
k3 = count_kmers(sequence, k)

k = 4
print("the frequency of the k-mers of length", k)
print(count_kmers(sequence, k))
k4 = count_kmers(sequence, k)

k = 5
print("the frequency of the k-mers of length", k)
print(count_kmers(sequence, k))
k5 = count_kmers(sequence, k)

################################################################################################

def find_max_kmers(kmer_dictionary):
    # compute the maximum frequency of a kmer
    max_frequency = max(kmer_dictionary.values())

    # a list to hold the k-mers with the maximum frequency
    max_kmers = []

    # iterate through the input dictionary to find all k-mers with the max frequency
    for kmer, count in kmer_dictionary.items():
        if count == max_frequency:
            max_kmers.append(kmer)

    return max_kmers

################################################################################################

k = 3
print("the most frequent k-mer of length", k)
max_kmers_k3 = find_max_kmers(k3)
print("K-mers of length", k, "with the highest frequency:", max_kmers_k3)

k = 4
print("the most frequent k-mer of length", k)
max_kmers_k4 = find_max_kmers(k4)
print("K-mers of length", k, "with the highest frequency:", max_kmers_k4)

k = 5
print("the most frequent k-mer of length", k)
max_kmers_k5 = find_max_kmers(k5)
print("K-mers of length", k, "with the highest frequency:", max_kmers_k5)

################################################################################################
################################################################################################
################################################################################################
################################################################################################ PART 3.b
# 3. Given the following DNA string: TGCAGGACGTATGGCCGTTAGCTAACGGCCATAGATGCTATG
#
#    a. Write a function which counts kmers for any given length (k). What is the most frequent 3, 4 and 5nt length kmer found in the sequence above?
#    b. Write a function which identifies the longest kmer in the string above that is also a palindrome.

def find_max_kmers_palindrome(kmer_dictionary):
    # compute the maximum frequency of a kmer
    max_frequency = max(kmer_dictionary.values())

    # a list to hold the k-mers with the maximum frequency
    max_kmers = []

    # iterate through the input dictionary to find all k-mers with the max frequency
    for kmer, count in kmer_dictionary.items():
        if count == max_frequency and palindrome_sequence(kmer):
            max_kmers.append(kmer)

    return max_kmers

###############################################################################################
###############################################################################################

k = 3
max_kmers_k3 = find_max_kmers(k3)
print("the most frequent K-mers of length", k, "that are palindromes:", find_max_kmers_palindrome(k3))

k = 4
max_kmers_k4 = find_max_kmers(k4)
print("the most frequent K-mers of length", k, "that are palindromes:", find_max_kmers_palindrome(k4))

k = 5
max_kmers_k5 = find_max_kmers(k5)
print("the most frequent K-mers of length", k, "that are palindromes:", find_max_kmers_palindrome(k5))

###############################################################################################
###############################################################################################
###############################################################################################
############################################################################################### PART 4

# Create a class to wrap together all the functions used to answer questions 1-3. Demonstrate how you would use this class to answer all the questions above.
# Write a main script which utilizes the palindrome finder class. Take as an input any DNA string and return the longest palindrome kmer identified from that string


class Palindrome_Finder2:
    """A class for finding the most frequent palindromes in a DNA sequence."""

    def reverse_complement2(self, seq):
        """This function returns a reverse complement of a DNA sequence"""
        # Complement the DNA strand
        seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").replace("N", "n").replace("-", "-")
        seq = seq.upper()

        # Reverse the complemented DNA strand
        seq = seq[::-1]
        return seq

    # returns wither True or False
    def palindrome_sequence(self, seq):
        """This function checks if the given DNA sequence is a palindrome"""

        # Get the reverse complement of the sequence
        rev_comp = self.reverse_complement2(seq)

        # Check if the sequence is the same as its reverse complement
        return seq.upper() == rev_comp.upper()

    def count_kmers(self, sequence, k):
        """This function counts the occurrences of k-mers in a DNA sequence."""
        kmer_counts = {}
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1
        return kmer_counts

    def find_max_kmers(self, kmer_dictionary):
        """This function finds the k-mers having the highest frequency."""
        # compute the maximum frequency of a kmer
        max_frequency = max(kmer_dictionary.values())

        # a list to hold the k-mers with the maximum frequency
        max_kmers = []

        # iterate through the input dictionary to find all k-mers with the max frequency
        for kmer, count in kmer_dictionary.items():
            if count == max_frequency:
                max_kmers.append(kmer)

        return max_kmers

    def find_max_kmers_palindrome(self, kmer_dictionary):
        """This function finds the k-mers having the highest frequency, that are palindromes"""
        # compute the maximum frequency of a kmer
        max_frequency = max(kmer_dictionary.values())

        # a list to hold the k-mers with the maximum frequency
        max_kmers = []

        # iterate through the input dictionary to find all k-mers with the max frequency
        for kmer, count in kmer_dictionary.items():
            if count == max_frequency and self.palindrome_sequence(kmer):
                max_kmers.append(kmer)

        return max_kmers

###############################################################################################
###############################################################################################
# running the script with the class Palindrome_Finder2() to replicate the results that were obtained above

print("")
print("")
print("")
sequence = "TGCAGGACGTATGGCCGTTAGCTAACGGCCATAGATGCTATG"
print("To replicate the results obtained on the sequence", sequence, "with the class Palindrome_Finder2() :")

ANALYSIS = Palindrome_Finder2()

###############################################################################################
############################################################################################### if k-mer length is 3
k = 3
print("")
print("the results of the analysis, considering the length of a k-mer", k)
kmer_counts = ANALYSIS.count_kmers(sequence, k)
max_kmers = ANALYSIS.find_max_kmers(kmer_counts)
max_kmers_palindrome = ANALYSIS.find_max_kmers_palindrome(kmer_counts)
print("The most frequent k-mer of length ", k)
print(max_kmers)
print("The most frequent k-mer of length ", k, "that is palindrome")
print(max_kmers_palindrome)

###############################################################################################
############################################################################################### if k-mer length is 3
k = 4
print("")
print("the results of the analysis, considering the length of a k-mer", k)
kmer_counts = ANALYSIS.count_kmers(sequence, k)
max_kmers = ANALYSIS.find_max_kmers(kmer_counts)
max_kmers_palindrome = ANALYSIS.find_max_kmers_palindrome(kmer_counts)
print("The most frequent k-mer of length ", k)
print(max_kmers)
print("The most frequent k-mer of length ", k, "that is palindrome")
print(max_kmers_palindrome)

###############################################################################################
############################################################################################### if k-mer length is 3
k = 5
print("")
print("the results of the analysis, considering the length of a k-mer", k)
kmer_counts = ANALYSIS.count_kmers(sequence, k)
max_kmers = ANALYSIS.find_max_kmers(kmer_counts)
max_kmers_palindrome = ANALYSIS.find_max_kmers_palindrome(kmer_counts)
print("The most frequent k-mer of length ", k)
print(max_kmers)
print("The most frequent k-mer of length ", k, "that is palindrome")
print(max_kmers_palindrome)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
print("")
print("")
print("")
############################################################################################### shall we use the command line and
############################################################################################### shall we provide a sequence
############################################################################################### shall we provide a k
###############################################################################################
###############################################################################################

def main():
    parser = argparse.ArgumentParser(description="To analyze a DNA sequences")
    parser.add_argument("sequence", type=str, help="The DNA sequence to analyze")
    parser.add_argument("k", type=int, help="The length of k-mers to analyze")

    args = parser.parse_args()

    print("")
    print("")
    print("")

    sequence = args.sequence
    k = args.k

    print("the input sequence is : ")
    print(sequence)
    print("the input k is : ")
    print(k)

    # Initialize the class
    ANALYSIS2 = Palindrome_Finder2()

    # Perform operations
    print("Reverse Complement:", ANALYSIS2.reverse_complement2(sequence))
    print("Is palindrome:", ANALYSIS2.palindrome_sequence(sequence))

    kmer_counts = ANALYSIS2.count_kmers(sequence, k)
    print("")
    print("")
    print("")
    print("Kmer counts for the sequence :", sequence, "and a k-mer :", k)
    print(kmer_counts)
    print("The kmers with the highest frequency:", ANALYSIS2.find_max_kmers(kmer_counts))
    print("The kmers with the highest frequency that are palindromic:", ANALYSIS2.find_max_kmers_palindrome(kmer_counts))

###############################################################################################
############################################################################################### we print the longest palindrome
############################################################################################### irrespective of a particular k value
###############################################################################################

    def find_longest_palindrome_kmer(sequence):
        """Finds the longest palindromic k-mer in the sequence, irrespective of k value."""
        longest_palindrome = ""
        sequence_length = len(sequence)
        for k in range(1, sequence_length + 1):
            for i in range(sequence_length - k + 1):
                kmer = sequence[i:i+k]
                if palindrome_sequence(kmer):
                    longest_palindrome = kmer
        return longest_palindrome

    print("")
    print("")
    print("")
    print("The following analysis considers any k value :")
    print("The longest Palindromic Kmer, irrespective of a particular k-value is : ", find_longest_palindrome_kmer(sequence))

if __name__ == "__main__":
    main()

###############################################################################################
###############################################################################################
###############################################################################################
