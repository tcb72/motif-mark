import re
from operator import itemgetter
from itertools import *

class Gene:

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

    def get_gene_name(self):
        return(self.header.split(' ')[0][1:])

    def get_intron_exon_locations(self):
        exon_indexes = []
        intron_indexes = []
        for index,char in enumerate(self.sequence):
            #deal w/ beginning and end of sequence
            if char.islower():
                intron_indexes.append(index)
            else:
                exon_indexes.append(index)

        intron_groups=[]
        exon_groups = []
        for k, g in groupby(enumerate(intron_indexes), lambda x: x[0]-x[1]):
            intron_groups.append(list(map(itemgetter(1), g)))

        for k, g in groupby(enumerate(exon_indexes), lambda x: x[0]-x[1]):
            exon_groups.append(list(map(itemgetter(1), g)))


        intron_start_end = [(item[0],item[-1]) for item in intron_groups]
        exon_start_end = [(item[0],item[-1]) for item in exon_groups]


        return intron_start_end, exon_start_end

    def find_motif_locations(self, motifs):
        motif_count = {}
        motif_pos_info = {}
        # for each motif
        for motif in motifs:
            curr_motif_start_positions = []
            curr_motif_seqs = [i.lower() for i in motifs[motif]]
            for i in range(len(self.sequence)):
                curr_kmer = self.sequence[i:i+len(motif)]
                if curr_kmer.lower() in curr_motif_seqs:
                    curr_motif_start_positions.append(i)
                if len(curr_kmer) < len(motif):
                    break
            motif_pos_info[motif] = curr_motif_start_positions

        return(motif_pos_info)




    def get_gene_length(self):
        return(len(self.sequence))
