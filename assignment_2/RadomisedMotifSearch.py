import sys
import math
import random
from collections import Counter


class FastAreader:
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''

        self.fname = fname
        self.fileH = None

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as self.fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = self.fileH.readline()
            while not line.startswith('>'):
                line = self.fileH.readline()
            header = line[1:].rstrip()

            for line in self.fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header, sequence

class CommandLine():
    """
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    available within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    """

    def __init__(self, inOpts=None):
        """
        CommandLine constructor.
        There are three arguments that are passed to the class Consensus():
            -i iterations (int)
            -k motif length (int)
            -p pseudocount (float)
        """
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='python randomizedMotifSearch.py -i int --maxMotif int --cutoff int < input.fa > output.out'
        )
        # Be sure to go over the argument information again
        self.parser.add_argument('-i', type=int, action='store',
                                 help='number of iterations (int)')
        self.parser.add_argument('-k', type=int, action='store',
                                 help='motif length (int)')
        self.parser.add_argument('-p', type=int, action='store',
                                 help='pseudocount (float)')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class Genome:
    """
    This class helps in initializing and storing different attributes that can be used across several functions across the class. 
    
    """
    def __init__(self, i, k, p, input_sequences):
        """

        """
        self.iterations = i
        self.pseudo_counts = p
        self.kmer_size = k
        self.input_sequences = input_sequences
        self.best_score = 0
        self.best_motif = []
        self.null_distribution = self.null_model_distribution()
        
    def newMotif(self, profile):
        """
        Function: To get new set of motifs from the existing profiles.
        input: profile of a set of existing selected motifs. 
        
        return: new set of motifs that generated using the counts of bases from the profile. 
        """


        new_motif = []
        for seq in self.input_sequences:
            best_prod_score = 0
            best_motif = []

            for i in range(0, len(seq)-self.kmer_size+1):
                motif = seq[i:i+self.kmer_size]
                product_score = 1
                for pos, base in enumerate(motif):
                    count_from_profile = profile[base][pos]
                    product_score = product_score * count_from_profile
                if product_score > best_prod_score:
                    best_prod_score = product_score
                    best_motif = motif
                
            new_motif.append(best_motif)
        
        return new_motif

    def selectRandomKmer(self, input_sequence_list):
        """
        Function: select random kmers from all the reads in the fast file. 
        input: input sequence list from the fasta file.
        
        return: set of random motifs of size k.
        """
        
        rand_kmer_array = []
        
        for sequence in input_sequence_list:
            
            rand_index = random.randint(0, len(sequence)-self.kmer_size)
            
            random_kmer = sequence[rand_index:rand_index+self.kmer_size]
            rand_kmer_array.append(random_kmer)
        return rand_kmer_array
    
    def null_model_distribution(self):
        """
        Function that returns a dictionary of number of counts of every bases {A,C,G,T} in the complete genome sequence.
        """
        input_seq = self.input_sequences
        joined_input_seq = "".join(input_seq)
        col_counter = dict(Counter(joined_input_seq))

        null_distribution = {k: (v+self.pseudo_counts) / (len(joined_input_seq)+ self.pseudo_counts*4) for k, v in col_counter.items()}
        
        return null_distribution


    def make_profile_from_kmers(self, input_motifs):
        """
        The function calculates the distribution of bases for the provided input set of motifs. 

        input: set of motifs of kmer-size k. 
        return: a maxtrix of size 4xk. where the four rows are the bases {A,C,G,T} and the columns of size k each having the counts of bases.
        """
        profile_dict = {
            'A': [],
            'G': [],
            'C': [],
            'T': []
        }

        for i in range(self.kmer_size):
            col = []
            for j in range(len(input_motifs)):
                col.append(input_motifs[j][i])
            col_counter = dict(Counter(col))
            
            nset = len(self.input_sequences)
            profile_dict['A'].append(col_counter['A']+self.pseudo_counts if 'A' in col_counter else 0+self.pseudo_counts)
            profile_dict['C'].append(col_counter['C']+self.pseudo_counts if 'C' in col_counter else 0+self.pseudo_counts)
            profile_dict['T'].append(col_counter['T']+self.pseudo_counts if 'T' in col_counter else 0+self.pseudo_counts)
            profile_dict['G'].append(col_counter['G']+self.pseudo_counts if 'G' in col_counter else 0+self.pseudo_counts)
        return profile_dict
    
    def calcRelativeEntropy(self, input_motifs):
        """
        The function considers the null model and compares it with the experimental model inorder to calculate the relative entropy.
        The function also call the make_profile_from_kmers function to create profile from the provided input set of motifs. 

        Input: the set of motifs of kmer-size k. 

        return: the relative entropy score of the provided set of motifs and the profile of the motifs. 

        """
        relative_motif_score = 0

        profile = self.make_profile_from_kmers(input_motifs)
        for col in range(self.kmer_size):
            for base in self.null_distribution:
                if profile[base][col] != 0: 
                    pr = profile[base][col]/(len(self.input_sequences) + (4*self.pseudo_counts))
                    relative_motif_score += pr * math.log2(pr / self.null_distribution[base])
        return relative_motif_score, profile

    def randomMotifSearch(self, input_seqs):
        """
        The function iterates for self.iterations times and initially selects random motifs of kmers-size k, and creates profile to select the next set of motifs which will be compared with the 
        previous set of motifs using relative entropy score. The process is repeated until the best relative entropy score is found. 

        Input: input sequences of the whole fasta file. 
        return: prints the consensus and the best relative entropy score of the given input. 
        """
        
        for i in range(self.iterations):
            
            random_kmers = self.selectRandomKmer(input_seqs) # initialized kmers randomly with kmer-size -k
            best_score_motif, best_profile = self.calcRelativeEntropy(random_kmers) # calculate the score of the randomly initialized kmers. Consider it as the best_score. 
            selected_motifs = random_kmers
            while True:
                new_motifs = self.newMotif(best_profile) # select new motifs from the current profile
                new_score, new_profile = self.calcRelativeEntropy(new_motifs)
                
                if new_score > best_score_motif:
                    selected_motifs = new_motifs
                    best_score_motif = new_score
                    best_profile = new_profile

                    if best_score_motif > self.best_score:
                        self.best_score = best_score_motif
                        self.best_motif = selected_motifs
                else:
                    break

        print("Found Consencus {} with score - {}".format(self.printConsencus(self.best_motif), self.best_score))

    def printConsencus(self, input_motifs):
        """
        Function: Using the best motif found using randomsearch algorithm, the function returns the consensus from the profile of the motif. 

        Input: input_motif - is the best motif found from the randmom search algorithm.

        return: the consensus motif of kmer-size k. 
        """
        motif_set = input_motifs
        profile = self.make_profile_from_kmers(motif_set)
        concensus = ''
        for i in range(self.kmer_size):
            max_prob = 0
            max_base = 'A'
            for base, probs in profile.items():
                if probs[i] > max_prob:
                    max_prob = probs[i]
                    max_base = base
            concensus+= max_base
        return concensus

def main(inFile="", options = None):
    ''' Setup necessary objects, read data and print the final report.'''
    cl = CommandLine() # setup the command line
    sourceReader = FastAreader() # setup the Fasta reader Object
    sequenceList = []
    for head, seq in sourceReader.readFasta():       # reading the fast file. 
        sequenceList.append(seq)
    
    thisGenome = Genome(cl.args.i, cl.args.k, cl.args.p, sequenceList)
    thisGenome.randomMotifSearch(sequenceList)

if __name__ == "__main__":
    main()
    