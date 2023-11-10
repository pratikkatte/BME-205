import sys
import math

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
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse

        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )

        self.parser.add_argument('-l', '--minMotif', nargs='?', default=1, action='store',
                                 help='min kMer size ')
        self.parser.add_argument('-m', '--maxMotif', nargs='?', default=8, action='store',
                                 help='max kMer size ')
        self.parser.add_argument('-c', '--cutoff', nargs='?', type=float, default=.01, action='store',
                                 help='Zscore cutoff')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class Genome:
    """
    This class helps in storing different attributes that can be used across several functions across the class. 
    Functions:
        reverseComplement: to get the reverse complement of the input sequence. 

    """
    def __init__(self, minK, maxK, cutoff, input_sequence):
        """Construct the object to count and analyze kMer frequencies"""
        self.min_kmer = int(minK)
        self.max_kmer = int(maxK)
        self.motif_counts = dict() #dict to store motif_counts
        self.genome_sequence = input_sequence
        self.cutoff = float(cutoff)
        self.genome_size = len(self.genome_sequence)

    def reverseComplement(self, input_seq):
        """
        reverse complement of the input sequence. 

        :returns - returns the reverse complement of the input motif. 
        """

        rev_seq_str = input_seq[::-1]
        reversse_compliment_seq = ''
        for i in rev_seq_str:
            if i=='A':
                reversse_compliment_seq = reversse_compliment_seq + 'T'
            elif i=='T':
                reversse_compliment_seq = reversse_compliment_seq + 'A'
            elif i=='C':
                reversse_compliment_seq = reversse_compliment_seq + "G"
            else:
                reversse_compliment_seq = reversse_compliment_seq + "C"
        return reversse_compliment_seq

    def motifSearch(self):
        """
        Function to find all the k-mers for size 3 to 8. It also calculates the frequency of each kmer in the input sequence. 
        output: 
            motif_counts: stored as class object. a dictonary with key is the kmers and value are the counts of the kmers. Saved in class variable self.motifCounts. stores all the motifs and counts in self.motif_counts
        """
        seq = self.genome_sequence
        for j in range(1, self.max_kmer+1):
            for k in range(0, len(seq) - j+ 1):
                motif = seq[k:k + j]
                if motif in self.motif_counts.keys(): 
                    rc = self.reverseComplement(motif)
                    self.motif_counts[motif] += 1
                    self.motif_counts[rc] = self.motif_counts[motif]
                else:
                    self.motif_counts[motif] = 1
                    rc = self.reverseComplement(motif)
                    self.motif_counts[rc] = 1

    def calcProbab(self, e):
        """
        Calculates Probability of a kmer appearing in a sequence.

        input: expect score of kmer. 
        return: probability value. 
        """
        prob = e/self.genome_size
        return prob

    def zScore(self, s, n, p):
        """
        This function calculates the z-score for a motif.
        input: 
            s - counts of a motif. 
            n - genomic size
            p - probability. 

        :return - returns the z-score if standard deviation is not zero, else returns 1
        """
        sd = math.sqrt(n*p*(1-p))
        if sd != 0:
            return (s - (n * p)) / sd
        else:
            return 1

    def expectedScore(self, motif):
        """
        Calculates the expected score of motif using HMM(2). It is done using the derived equation that consists of three terms.
        1st term - Considers the motif excluding the last nucleotide. 
        2nd term - considers the motif excluding the first nucleotide. 
        3rd term (denominator) - considers the middle part of the motif, excluding the first and last nucleotide.

        return: expect score
        The function returns zero if the motif is not found in the motif_counts dictonary. 
        """
        if motif in self.motif_counts:
            first_part =  motif[:-1]
            second_part = motif[1:]
            denominator_part = motif[1:-1]

            if first_part in self.motif_counts and second_part in self.motif_counts and denominator_part in self.motif_counts:
                first_motif_counts = self.motif_counts[first_part]
                second_motif_counts = self.motif_counts[second_part]
                denominator_motif_counts = self.motif_counts[denominator_part]
                est_score = (first_motif_counts * second_motif_counts)/(denominator_motif_counts)
                return est_score
            else:
                return 0

    def motifScores(self):
        """
        Calls the function expectedScore(), zScore(), and calcProbab() to calculate expected_score and z-score for every motif that is stored in motif_counts dictionary. 
        
        :return - a dictionary (final_dict) which consists of all the motifs as 'key' and its 'values' are [reverse component, expected score and z-score].
        """
        is_motif_calculated = dict()
        final_dict = dict()
        for motif, counts in self.motif_counts.items():
            rc_motif = self.reverseComplement(motif)
            if motif not in is_motif_calculated and rc_motif not in is_motif_calculated:
                kmers = [motif, rc_motif]
                exp_score = self.expectedScore(motif)
                s_count = self.motif_counts[motif]
                prob = self.calcProbab(exp_score)
                z_score = self.zScore(s_count, self.genome_size, prob)
                if z_score <= self.cutoff:
                    kmers.sort()
                    final_dict[kmers[0]] = [kmers[1], counts, exp_score, z_score]
                is_motif_calculated[motif] = True
                is_motif_calculated[rc_motif] = True
        return final_dict
    
    def display_scores(self, input_dict):
        """
        Displays all the motifs with size in decending orders and within every motif-size, sorts the z-score in ascending order. 


        :input- a dictionary that consists of motifs as 'keys' and values as the following [reverse component, counts, expected score, and z-score]
        :return- prints the output in required format. 
        """
        sorted_dict = dict(sorted(input_dict.items(), key=lambda item: (-len(item[0]), item[1][-1])))

        print("N = {}".format(self.genome_size-8))
        for key,values in sorted_dict.items():
            print("{}:{}\t {}\t {}\t {}".format(key, values[0], values[1], values[2], values[3]))

def main(inFile="", options = None):

    ''' Setup necessary objects, read data and print the final report.'''
    cl = CommandLine() # setup the command line
    sourceReader = FastAreader() # setup the Fasta reader Object
    sequenceList = []
    for head, seq in sourceReader.readFasta():       # reading the fast file. 
        sequenceList.extend(seq)
    
    sequenceList = "".join(sequenceList) # concatenating all the sequencing into one single string. 

    thisGenome = Genome(cl.args.minMotif, cl.args.maxMotif, cl.args.cutoff, sequenceList) # setup a Genome object

    k = thisGenome.motifSearch()
    output = thisGenome.motifScores()
    thisGenome.display_scores(output)

if __name__ == "__main__":  
    main()

