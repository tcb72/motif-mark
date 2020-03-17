import cairocffi as cairo
from gene import Gene
from itertools import product
import random
import seaborn as sns
import argparse

def multiline_to_dict_fasta(file):
    fasta_dict = dict()
    with open(file) as f:
        for index,line in enumerate(f):
            if (line.startswith('>')):
                curr_header = line.strip()
                fasta_dict[curr_header] = ''
            else:
                fasta_dict[curr_header] += line.strip()
    return(fasta_dict)

def get_motifs(file):
    with open(file) as f:
        raw_motifs = [i.strip() for i in f.readlines()]

    ambiguous_bases = {'y':['c','t'], 'u':['u','t']}
    motif_dict = {}
    for motif in raw_motifs:
        motif_char_list = []
        motif_chars = list(motif)
        for char in motif_chars:
            if char.lower() in ambiguous_bases:
                motif_char_list.append(ambiguous_bases[char.lower()])
            else:
                motif_char_list.append([char])
        if motif.isupper():
            curr_motifs = [''.join(i).upper() for i in list(product(*motif_char_list))]
        else:
            curr_motifs = [''.join(i) for i in list(product(*motif_char_list))]
        motif_dict[motif] = curr_motifs
    return(motif_dict)

def draw_surface(fasta_dict, num_motifs):
    num_genes = len(fasta_dict)
    # width is length of longest sequence plus 100 (for some extra white space)
    WIDTH = len(sorted(fasta_dict.values(), key=len)[-1]) + 100
    # height is the y_offset * number of genes
    HEIGHT = 100*num_genes + 20*num_motifs
    surface = cairo.SVGSurface("plot.svg", WIDTH, HEIGHT)
    context = cairo.Context(surface)
    return(context)

def draw_introns(introns, X_OFFSET, y_offset, Y_INITIAL):
    for intron in introns:
        intron_start = intron[0]
        intron_end = intron[1]
        context.set_source_rgb(0, 0, 0)
        context.move_to(X_OFFSET+intron_start,Y_INITIAL+y_offset)
        context.line_to(X_OFFSET+intron_end,Y_INITIAL+y_offset)
        context.stroke()

def draw_exons(exons, X_OFFSET, y_offset, Y_INITIAL, RECTANGLE_HEIGHT):
    for exon in exons:
        exon_start = exon[0]
        exon_end = exon[1]
        context.set_source_rgb(0, 0, 0)
        context.rectangle(X_OFFSET+exon_start, Y_INITIAL-(RECTANGLE_HEIGHT/2)+y_offset, exon_end-exon_start, RECTANGLE_HEIGHT)
        context.stroke()

def draw_motifs(palette, motif_locations,X_OFFSET, y_offset, Y_INITIAL):
    for index,motif in enumerate(motif_locations):
        context.set_source_rgb(*palette[index])
        for pos in motif_locations[motif]:
            context.move_to(X_OFFSET + pos, Y_INITIAL + y_offset + 10)
            context.line_to(X_OFFSET + pos, Y_INITIAL + y_offset - 10)
            context.stroke()

def draw_legend(palette, X_OFFSET, y_offset,motifs):
    legend_offset = 0
    for index,rgb in enumerate(palette):
        context.move_to(X_OFFSET+15,y_offset+legend_offset+3)
        context.set_source_rgb(0, 0, 0)
        context.show_text(motifs[index])
        context.set_source_rgb(rgb[0],rgb[1],rgb[2])
        context.move_to(X_OFFSET,y_offset+legend_offset)
        context.line_to(X_OFFSET+10, y_offset+legend_offset)
        context.stroke()
        legend_offset += 20



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta')
    parser.add_argument('--motifs')
    args = parser.parse_args()

    #'/home/tcb/Downloads/Figure_1.fasta'
    #'/home/tcb/Downloads/Fig_1_motifs.txt'
    fasta = multiline_to_dict_fasta(args.fasta)
    motifs = get_motifs(args.motifs)


    X_OFFSET = 75
    Y_INITIAL = 25
    RECTANGLE_HEIGHT = 30
    palette = sns.color_palette(None, len(motifs))
    context = draw_surface(fasta, len(motifs))
    context.set_line_width(1)
    context.save()
    context.set_source_rgb(1, 1, 1)
    context.paint()
    context.restore()

    y_offset = 0
    for header,sequence in fasta.items():
        gene = Gene(header,sequence)
        gene_name = gene.get_gene_name()
        gene_length = gene.get_gene_length()
        introns, exons = gene.get_intron_exon_locations()
        motif_locations = gene.find_motif_locations(motifs)
        context.set_source_rgb(0, 0, 0)
        context.move_to(20,Y_INITIAL+y_offset+5)
        context.show_text(gene_name)

        draw_introns(introns,X_OFFSET, y_offset, Y_INITIAL)
        draw_exons(exons, X_OFFSET, y_offset, Y_INITIAL, RECTANGLE_HEIGHT)
        draw_motifs(palette,motif_locations, X_OFFSET, y_offset, Y_INITIAL)

        y_offset += 100

    draw_legend(palette,X_OFFSET,y_offset,list(motifs.keys()))
