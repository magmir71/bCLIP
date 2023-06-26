#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Author : Aleksei Mironov
# Company: Mihaela Zavolan, Biozentrum, Basel
# This script is part of the Zavolan lab ZARP pipeline.
# -----------------------------------------------------------------------------

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd

def main():
    """ Get genomic segmentation out of the .gtf file"""

    __doc__ = "Get genomic segmentation out of the .gtf file"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_gtf",
                        dest="input_gtf_file_path",
                        help="Path to the input gtf annotation file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--input_genome_fai",
                        dest="input_genome_fai_file_path",
                        help="Path to the input genome fai index file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--output_exon_intron_gs",
                        dest="output_exon_intron_gs_file_path",
                        help="path for the output file output_exon_intron_genome_segmentation.bed",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--output_first_middle_last_gs",
                        dest="output_first_middle_last_gs_file_path",
                        help="path for the output file output_first_middle_last_genome_segmentation.bed",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--temp_dir",
                    dest="temp_dir_path",
                    help="path to the temporary directory",
                    required=True,
                    metavar="FILE")
    parser.add_argument("--transcript_ends_length",
                    dest="transcript_ends_length",
                    help="FIRST and LAST sub-segments have the specified length. Used for the analysis of 5'-bias and 3'-bias in expression",
                    required=True,
                    metavar="FILE")
    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    temp_dir_path = options.temp_dir_path+'/'
    
    # create output directories if they don't exist
    for file_path in [options.output_exon_intron_gs_file_path, options.output_first_middle_last_gs_file_path, temp_dir_path]:
        outdir = '/'.join((file_path.split('/')[:-1]))+'/'
        out = subprocess.check_output('mkdir -p '+outdir, shell=True)
    
    input_gtf = pd.read_csv(options.input_gtf_file_path,delimiter="\t",index_col=None,header=None)
    
    for sel_strand in ['+','-']:

        cur_input_gtf = input_gtf.loc[input_gtf[6]==sel_strand].reset_index(drop=True).copy()

        genes = cur_input_gtf.loc[cur_input_gtf[2]=='gene'].reset_index(drop=True)
        genes['gene_id'] = genes[8].str.split('gene_id "',expand=True)[1].str.split('";',expand=True)[0]
        genes['score'] = 0
        genes = genes[[0,3,4,'gene_id','score',6]]
        genes[3] = genes[3]-1
        genes = genes.sort_values([0,3,4],ascending=True).reset_index(drop=True)
        genes.to_csv(temp_dir_path+'genes_original.bed', sep=str('\t'),header=False,index=None)

        # putting overlapping genes into disjoint regions
        command = 'bedops --partition '+temp_dir_path+'genes_original.bed > '+temp_dir_path+'genes_partion.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools intersect -loj -sorted -a '+temp_dir_path+'genes_partion.bed -b '+temp_dir_path+\
        'genes_original.bed | bedtools groupby -g 1,2,3,8,9 -c 7 -o distinct > '+temp_dir_path+'genes_disjoint.bed'
        out = subprocess.check_output(command, shell=True)
        genes_disjoint = pd.read_csv(temp_dir_path+'genes_disjoint.bed',delimiter="\t",index_col=None,header=None)
        genes_disjoint[[0,1,2,5,3,4]].to_csv(temp_dir_path+'genes_disjoint.bed', sep=str('\t'),header=False,index=None)

        exons = cur_input_gtf.loc[cur_input_gtf[2]=='exon'].reset_index(drop=True)
        exons['exon_id'] = exons[8].str.split('exon_id "',expand=True)[1].str.split('";',expand=True)[0]
        exons['score'] = 0
        exons = exons[[0,3,4,'exon_id','score',6]]
        exons[3] = exons[3]-1
        exons = exons.sort_values([0,3,4],ascending=True).reset_index(drop=True)
        exons.to_csv(temp_dir_path+'exons.bed', sep=str('\t'),header=False,index=None)

        # get introns
        exons = cur_input_gtf.loc[cur_input_gtf[2]=='exon'].reset_index(drop=True)
        exons['transcript_id'] = exons[8].str.split('transcript_id "',expand=True)[1].str.split('";',expand=True)[0]
        exons = exons[[0,3,4,'transcript_id',6]]
        exons = exons.sort_values([0,6,'transcript_id',3,4],ascending=True).reset_index(drop=True)
        exons = exons.values
        introns = []
        intron_start, intron_end, prev_transcript = None, None, ''

        for exon in exons:
            if (intron_start is None) or (prev_transcript!=exon[3]):
                intron_start = exon[2]
                prev_transcript = exon[3]
            else:
                intron_end = exon[1]
                introns.append([exon[0],intron_start,intron_end,exon[3],0,exon[4]])
                intron_start = exon[2]
                prev_transcript = exon[3]
        introns = pd.DataFrame(introns)
        introns[2] = introns[2]-1
        introns = introns.sort_values([0,1,2],ascending=True).reset_index(drop=True)
        introns.to_csv(temp_dir_path+'introns.bed', sep=str('\t'),header=False,index=None)

        genome_fai = pd.read_csv(options.input_genome_fai_file_path,delimiter="\t",index_col=None,header=None)
        genome_fai['start'] = 0
        genome_fai['score'] = 0
        genome_fai['chr_name'] = genome_fai[0]
        genome_fai['strand'] = sel_strand
        genome_fai[[0,'start',1,'chr_name','score','strand']].to_csv(temp_dir_path+'chromosomes.bed', sep=str('\t'),header=False,index=None)

        command = 'bedtools merge -i '+temp_dir_path+'exons.bed > '+temp_dir_path+'exons_collapsed.bed'
        out = subprocess.check_output(command, shell=True)

        command = 'bedtools merge -i '+temp_dir_path+'introns.bed > '+temp_dir_path+'introns_collapsed.bed'
        out = subprocess.check_output(command, shell=True)

        command = 'bedtools subtract -a '+temp_dir_path+'genes_disjoint.bed -b '+temp_dir_path+'introns_collapsed.bed > '+temp_dir_path+'alw_exonic.bed'
        out = subprocess.check_output(command, shell=True)

        command = 'bedtools subtract -a '+temp_dir_path+'genes_disjoint.bed -b '+temp_dir_path+'exons_collapsed.bed > '+temp_dir_path+'alw_intronic.bed'
        out = subprocess.check_output(command, shell=True)

        command = 'bedtools subtract -a '+temp_dir_path+'genes_disjoint.bed -b '+temp_dir_path+'alw_exonic.bed > '+temp_dir_path+'temp.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools subtract -a '+temp_dir_path+'temp.bed -b '+temp_dir_path+'alw_intronic.bed > '+temp_dir_path+'exonic_or_intronic.bed'
        out = subprocess.check_output(command, shell=True)


        command = 'bedtools merge -i '+temp_dir_path+'genes_disjoint.bed > '+temp_dir_path+'genes_collapsed.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools subtract -a '+temp_dir_path+'chromosomes.bed -b '+temp_dir_path+'genes_collapsed.bed > '+temp_dir_path+'intergenic.bed'
        out = subprocess.check_output(command, shell=True)

        exonic_or_intronic = pd.read_csv(temp_dir_path+'exonic_or_intronic.bed',delimiter="\t",index_col=None,header=None)
        alw_exonic = pd.read_csv(temp_dir_path+'alw_exonic.bed',delimiter="\t",index_col=None,header=None)
        alw_intronic = pd.read_csv(temp_dir_path+'alw_intronic.bed',delimiter="\t",index_col=None,header=None)
        intergenic = pd.read_csv(temp_dir_path+'intergenic.bed',delimiter="\t",index_col=None,header=None)

        alw_exonic['segment'] = 'alw_exonic'
        alw_intronic['segment'] = 'alw_intronic'
        exonic_or_intronic['segment'] = 'exonic_or_intronic'
        intergenic['segment'] = 'intergenic'

        gene_segments = pd.concat([alw_exonic,alw_intronic,exonic_or_intronic,intergenic])
        gene_segments = gene_segments.sort_values([0,1,2],ascending=True).reset_index(drop=True)
        gene_segments[3] = gene_segments['segment']+';'+gene_segments[3]
        gene_segments[3] = gene_segments.groupby([0,1,2,4,5])[3].transform(lambda x: ' OR '.join(x))
        gene_segments = gene_segments[[0,1,2,3,4,5]].drop_duplicates()

        gene_segments[[0,1,2,3,4,5]].to_csv(temp_dir_path+'genome_segmentation_'+sel_strand+'.bed', sep=str('\t'),header=False,index=None)

    genome_segmentation_plus = pd.read_csv(temp_dir_path+'genome_segmentation_+.bed',delimiter="\t",index_col=None,header=None)
    genome_segmentation_minus = pd.read_csv(temp_dir_path+'genome_segmentation_-.bed',delimiter="\t",index_col=None,header=None)

    # Now, to get antisense regions, we overlap genic regions on one strand with intergenic on the other
    genome_segmentation_plus.loc[~genome_segmentation_plus[3].str.contains('intergenic')].to_csv(temp_dir_path+'genic_segments_plus.bed', sep=str('\t'),header=False,index=None)
    genome_segmentation_minus.loc[genome_segmentation_minus[3].str.contains('intergenic')].to_csv(temp_dir_path+'intergenic_segments_minus.bed', sep=str('\t'),header=False,index=None)
    command = 'bedtools intersect -wb -a '+temp_dir_path+'intergenic_segments_minus.bed -b '+temp_dir_path+'genic_segments_plus.bed -S > '+temp_dir_path+'antisense_segments_minus.bed'
    out = subprocess.check_output(command, shell=True)

    genome_segmentation_minus.loc[~genome_segmentation_minus[3].str.contains('intergenic')].to_csv(temp_dir_path+'genic_segments_minus.bed', sep=str('\t'),header=False,index=None)
    genome_segmentation_plus.loc[genome_segmentation_plus[3].str.contains('intergenic')].to_csv(temp_dir_path+'intergenic_segments_plus.bed', sep=str('\t'),header=False,index=None)
    command = 'bedtools intersect -wb -a '+temp_dir_path+'intergenic_segments_plus.bed -b '+temp_dir_path+'genic_segments_minus.bed -S > '+temp_dir_path+'antisense_segments_plus.bed'
    out = subprocess.check_output(command, shell=True)

    # now substract antisense regions from intergenic regions
    command = 'bedtools subtract -a '+temp_dir_path+'intergenic_segments_plus.bed -b '+temp_dir_path+'antisense_segments_plus.bed > '+temp_dir_path+'refined_intergenic_segments_plus.bed'
    out = subprocess.check_output(command, shell=True)
    command = 'bedtools subtract -a '+temp_dir_path+'intergenic_segments_minus.bed -b '+temp_dir_path+'antisense_segments_minus.bed > '+temp_dir_path+'refined_intergenic_segments_minus.bed'
    out = subprocess.check_output(command, shell=True)

    antisense_segments_minus = pd.read_csv(temp_dir_path+'antisense_segments_minus.bed',delimiter="\t",index_col=None,header=None)
    antisense_segments_minus[9] = 'ANTISENSE:'+antisense_segments_minus[9]
    antisense_segments_minus = antisense_segments_minus[[0,1,2,9,4,5]]
    antisense_segments_minus.columns = [0,1,2,3,4,5]
    antisense_segments_plus = pd.read_csv(temp_dir_path+'antisense_segments_plus.bed',delimiter="\t",index_col=None,header=None)
    antisense_segments_plus[9] = 'ANTISENSE:'+antisense_segments_plus[9]
    antisense_segments_plus = antisense_segments_plus[[0,1,2,9,4,5]]
    antisense_segments_plus.columns = [0,1,2,3,4,5]

    refined_intergenic_segments_minus = pd.read_csv(temp_dir_path+'refined_intergenic_segments_minus.bed',delimiter="\t",index_col=None,header=None)
    refined_intergenic_segments_plus = pd.read_csv(temp_dir_path+'refined_intergenic_segments_plus.bed',delimiter="\t",index_col=None,header=None)
    genic_segments_minus = genome_segmentation_minus.loc[~genome_segmentation_minus[3].str.contains('intergenic')].reset_index(drop=True)
    genic_segments_plus = genome_segmentation_plus.loc[~genome_segmentation_plus[3].str.contains('intergenic')].reset_index(drop=True)

    genome_segmentation = pd.concat([genic_segments_plus,refined_intergenic_segments_plus,antisense_segments_plus,
                                     genic_segments_minus,refined_intergenic_segments_minus,antisense_segments_minus]).sort_values([0,1,2],ascending=True).reset_index(drop=True)
    genome_segmentation.to_csv(options.output_exon_intron_gs_file_path, sep=str('\t'),header=False,index=None)
    
    ### NOW, we get the sub-segmentation of alw_exonic segments into three sub-segments - FIRST, MIDDLE, LAST, to analyze 5'-bias and 3'-bias in expression. 
    ### We take segments longer than 3*<transcript_ends_length> nt, separate first <transcript_ends_length> nts, last <transcript_ends_length> nts, and the rest
    transcript_ends_length = int(options.transcript_ends_length)
    genome_segmentation['len'] = genome_segmentation[2]-genome_segmentation[1]
    gr = genome_segmentation.loc[genome_segmentation[3].str.contains('alw_exonic')].groupby([3]).agg({'len':sum}).reset_index().rename(columns={'len':'total_len'})
    sel_genome_segmentation = pd.merge(genome_segmentation,gr.loc[gr['total_len']>3*transcript_ends_length],how='inner',on=[3]).sort_values([0,5,3]).reset_index(drop=True)

    thr = transcript_ends_length
    a = []
    prev_segment = None
    prev_prefix = None
    prev_pos = 0
    for elem in sel_genome_segmentation.values:
        cur_segment = elem[3]
        strand = elem[5]
        step1, step2, step3 = 0, 0, 0
        if cur_segment!=prev_segment:
            step1 = min(thr,elem[6])
            prefix = 'FIRST>' if strand == '+' else 'LAST>'
            l = list(elem[[0,1]])+[elem[1]+step1,prefix+elem[3]]+list(elem[[4,5]])
            a.append(l)
            # if there is still a left part of the region, save it as part of the middle region 
            if step1<elem[6]:
                prefix = 'MIDDLE>'
                step2 = min(elem[6]-step1,elem[7]-2*thr)
                l = [elem[0],elem[1]+step1,elem[1]+step1+step2,prefix+elem[3]]+list(elem[[4,5]])
                a.append(l)
                # if there is still a left part of the region, save it as the last/first
                if step1+step2<elem[6]:
                    prefix = 'LAST>' if strand == '+' else 'FIRST>'
                    step3 = min(elem[6]-step1-step2,thr)
                    l = [elem[0],elem[1]+step1+step2,elem[1]+step1+step2+step3,prefix+elem[3]]+list(elem[[4,5]])
                    a.append(l)
            prev_pos = step1+step2+step3
        else:
            if prev_prefix in ['FIRST>','LAST>']:
                step1 = min(thr-prev_pos,elem[6]) if prev_pos<thr else min(elem[7]-prev_pos,elem[6])
                if step1>0:
                    prefix = 'FIRST>' if prev_prefix == 'FIRST>' else 'LAST>'
                    l = list(elem[[0,1]])+[elem[1]+step1,prefix+elem[3]]+list(elem[[4,5]])
                    a.append(l)
                # if there is still a left part of the region, save it as part of the middle region
                if step1<elem[6]:
                    prefix = 'MIDDLE>'
                    step2 = min(elem[6]-step1,elem[7]-prev_pos-thr)
                    l = [elem[0],elem[1]+step1,elem[1]+step1+step2,prefix+elem[3]]+list(elem[[4,5]])
                    a.append(l)
                    # if there is still a left part of the region, save it as the last/first
                    if step1+step2<elem[6]:
                        prefix = 'LAST>' if strand == '+' else 'FIRST>'
                        step3 = min(elem[6]-step1-step2,thr)
                        l = [elem[0],elem[1]+step1+step2,elem[1]+step1+step2+step3,prefix+elem[3]]+list(elem[[4,5]])
                        a.append(l)
            elif prev_prefix == 'MIDDLE>':
                step2 = min(elem[6]-step1,elem[7]-prev_pos-thr)
                if step2>0:
                    prefix = 'MIDDLE>'
                    l = list(elem[[0,1]])+[elem[1]+step1+step2,prefix+elem[3]]+list(elem[[4,5]])
                    a.append(l)
                if step2<elem[6]:
                    prefix = 'LAST>' if strand == '+' else 'FIRST>'
                    step3 = min(elem[6]-step2,thr)
                    l = [elem[0],elem[1]+step1+step2,elem[1]+step1+step2+step3,prefix+elem[3]]+list(elem[[4,5]])
                    a.append(l)
            prev_pos = prev_pos+step1+step2+step3
        prev_segment = cur_segment
        prev_prefix = prefix

    pd.DataFrame(a).sort_values([0,1,2],ascending=True).to_csv(options.output_first_middle_last_gs_file_path, sep=str('\t'),header=False,index=None)
    
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)