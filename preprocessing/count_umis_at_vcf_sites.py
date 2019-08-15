import pysam
import pandas as pd
import numpy as np
import argparse
import scipy.weave
import sys
from scipy.weave import converters

def get_umi_histogram(read_forward, read_reverse, match_options):
    umi_hist = {}
    for qname in read_forward:
        if read_reverse.has_key(qname):
            umi = qname.split(':')[-1]
            if umi.count('N') == 0:
                umi_hist[umi] = umi_hist.get(umi,0) + 1
    return umi_hist

def get_umi_from_qname(qname):
    seq_id = qname.split(':')
    umi = seq_id[-2] + seq_id[-1]
    return umi

def get_reads_for_target(target, bam, match_options, do_filter = False):
    """
    Returns dictionary of id -> bam_read for forward reads and reverse reads.
    Performs filtering using filter parameters in match_options
    """
    read_iter = bam.fetch(target[0], target[1], target[2])
    n = 0
    nall = 0
    forward_reads = {}
    reverse_reads = {}
    for read in read_iter:
        nall += 1
        if do_filter:
            if filter_read(read, target, match_options):
                continue
        if read.is_read1:
            forward_reads[read.qname] = read
        elif read.is_read2:
            reverse_reads[read.qname] = read
        else:
            1/0

        n += 1
    return forward_reads, reverse_reads

def get_majority_vote(nuc_dict = {}):
    m = -1
    basecall = 'nocall'
    for nuc, count in nuc_dict.iteritems():
        if count > m:
            m = count
            basecall = nuc
    return basecall

def count_umis(vcf_file, bamfile, match_options):
    bamfiles = [f for f in bamfile.split(',') if f not in [',','']]
    assert len(bamfiles) == 1,"Not output BAM filename to dataframe, that's why"
    bam = {}
    for bf in set(bamfiles):
        bam[bf] = pysam.AlignmentFile(bf, "rb")


    # output dataframe
    columns = ['CHROM',
              'POS',
              'ALT',
              'REF',
        #      'bamfile', 
              'cell_id',
              'A',
              'C',
              'G',
              'T',
              'N',
              'no_umi_A',
              'no_umi_C',
              'no_umi_G',
              'no_umi_T',
              'no_umi_N'
              ]
    dfs = []

    vcf = pysam.VariantFile(vcf_file,'r')

    count = 0

    for vcf_rec in vcf.fetch(region = match_options.vcf_region):
        if not (len(vcf_rec.alts[0])==1 and len(vcf_rec.ref)==1): # only consider SNPs 
            print "Skipping non-SNP",vcf_rec.chrom,vcf_rec.pos,vcf_rec.ref,",".join(vcf_rec.alts)
            continue
        #if len(vcf_rec.alts)>1: # skip multiallelic variants
        #    print "Skipping multiallelic variants %s:%d" % (vcf_rec.chrom, vcf_rec.pos)
        #    continue
 

        # target is a window around the target ASE variant such that also the mate read is fetched.
        # what if the variant 
        target = [ vcf_rec.chrom, max(vcf_rec.pos - match_options.window_around_ase_variant, 0), vcf_rec.pos + match_options.window_around_ase_variant]
        for bamfilename in bam:
            num_qual_skipped = 0
            num_no_cell_barcode = 0
            num_no_umi = 0
            qnames_seen = {}
            
            # cell -> molecule -> allele_count dictionary
            cell_to_molecule_to_allele_count = {}

    
            for pileupcolumn in bam[bamfilename].pileup(vcf_rec.chrom, vcf_rec.start, vcf_rec.stop,**{'max_depth':1000000} ): #FIXME
                
                if pileupcolumn.reference_pos == vcf_rec.start: #pileup position is zero-based #FIXME
                    # print ("\tcoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                    unmatched_read = 0
                    # print "\tpc:",pileupcolumn.pos
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            # query position is None if is_del or is_refskip is set.
                            # print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position])) 
                            if pileupread.alignment.mapq<match_options.min_mapping_qual:
                                continue
                            qqual = pileupread.alignment.query_qualities[pileupread.query_position]
                            if qqual < match_options.min_base_qual:
                                num_qual_skipped += 1
                                continue
                            qname = pileupread.alignment.query_name
                            qnames_seen[qname] = qnames_seen.get(qname, 0) + 1
                            assert qnames_seen[qname]  <= 2
                            qbase = pileupread.alignment.query_sequence[pileupread.query_position]

                            try:
                                cell_barcode = pileupread.alignment.get_tag('CB')
                            except KeyError:
                                num_no_cell_barcode += 1
                                continue
                            try:
                                umi = pileupread.alignment.get_tag('UB')
                            except KeyError:
                                num_no_umi += 1
                                continue
                            if not cell_barcode in cell_to_molecule_to_allele_count:
                                cell_to_molecule_to_allele_count[cell_barcode] = {}
                            molecule_to_allele_count = cell_to_molecule_to_allele_count[cell_barcode]
                            if not umi in molecule_to_allele_count:
                                molecule_to_allele_count[umi]  = {'A':0,'C':0,'G':0,'T':0,'N':0}
                            molecule_to_allele_count[umi][qbase] += 1
            # now output for each variant a table for allele counts in each cell
            if len(cell_to_molecule_to_allele_count) < match_options.min_number_of_cells_with_molecules:
                # print "Skipping because of number of cells",vcf_rec.chrom,vcf_rec.pos,vcf_rec.ref,",".join(vcf_rec.alts)
                continue
            cell_rows = []
            total_molecule_count = 0
            total_non_ref_molecule_count = 0
            for cell in cell_to_molecule_to_allele_count:
                rd = {'CHROM':vcf_rec.chrom,
                      'POS':vcf_rec.pos,
                      'ALT':",".join(vcf_rec.alts),
                      'REF':vcf_rec.ref,
                      #'bamfile' : bamfilename, 
                      'cell_id' : cell,
                      'A' : 0,
                      'C' : 0,
                      'G' : 0,
                      'T' : 0,
                      'N' : 0,
                      'no_umi_A' : 0,
                      'no_umi_C' : 0,
                      'no_umi_G' : 0,
                      'no_umi_T' : 0,
                      'no_umi_N' : 0
                      }
            
                
                pdict = cell_to_molecule_to_allele_count[cell]
                print(pdict)
                total_molecule_count += len(pdict)
                for umi in pdict:
                    basecall = get_majority_vote(pdict[umi])
                    rd[basecall] += 1
                    # this counts alleles without taking into account the umi's, for purpose of accuracy comparison.
                    for basecall in ['A','C','G','T','N']:
                        rd['no_umi_'+basecall] += pdict[umi][basecall]
                total_non_ref_molecule_count += rd['A']+rd['C']+rd['G']+rd['T'] - rd[vcf_rec.ref]
                cell_rows.append(rd)
                # print(rd)
            if total_molecule_count < match_options.min_total_molecule_count:
                continue
            if total_non_ref_molecule_count < match_options.min_total_nonref_molecule_count:
                print "Skipping non-variable SNP %s:%d:%s:%s" % (vcf_rec.chrom, vcf_rec.pos, vcf_rec.ref, ",".join(vcf_rec.alts))
                continue
            fraction_nonref = float(total_non_ref_molecule_count)/float(total_molecule_count)
            if fraction_nonref < match_options.min_fraction_nonref_molecule_count:
                print "Skipping SNP with fraction nonref =%1.2f < %1.2f %s:%d:%s:%s" % (fraction_nonref, match_options.min_fraction_nonref_molecule_count, vcf_rec.chrom, vcf_rec.pos, vcf_rec.ref, ",".join(vcf_rec.alts))
                continue
            # print(cell_rows)
            # exit()
            print(vcf_rec.alts)
            df = pd.DataFrame(cell_rows)[columns]
            # print(df)
            # exit()
            dfs.append(df)

            print vcf_rec.chrom,vcf_rec.pos,vcf_rec.ref,",".join(vcf_rec.alts),
            print "\tNumber of cells observed:",len(cell_to_molecule_to_allele_count),
            print "\tNumber of low-quality nucleotides:",num_qual_skipped,
            print "\tNumber of read-pairs seen in pileup:",len(qnames_seen.keys()),
            print "\tTotal molecule count",total_molecule_count,"non-ref molecule count",total_non_ref_molecule_count,
            print "\tA,C,G,T:",df['A'].sum(), df['C'].sum(), df['G'].sum(), df['T'].sum()
           
                # break
        #if len(dfs)>0:
        #    break
    print(len(dfs), " length dfs")
    if len(dfs)>0: 
        output_df = pd.concat(dfs, ignore_index = True)
        output_df.to_csv(match_options.output_file, sep = "\t", index=False)
    else:
        output_df = None
    
    for bamfilename in bam: 
        bam[bamfilename].close()
    
    return output_df

def main():
    #python 6_count_umis_at_vcf_sites.py --input 6samples_2.recode.chr1.vcf.gz --bamfile ../sortedBam.bam --vcf_region 1:100000-1100000 --output_file testOutput.txt
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf_file", help="mip design file", required = True)
    parser.add_argument("--bamfile", help="indexed BAM file", required = True)
    parser.add_argument("--output_file", help="tab-separated output file with molecule counts", required = True)

    parser.add_argument("--window_around_ase_variant", type = int, default = 1000)
    parser.add_argument("--min_base_qual", type = int, default = 15)
    parser.add_argument("--min_number_of_cells_with_molecules", type = int, default = 1)
    parser.add_argument("--min_total_molecule_count", type = int, default = 10)
    parser.add_argument("--min_total_nonref_molecule_count", type = int, default = 2)
    parser.add_argument("--min_fraction_nonref_molecule_count", type = float, default = 0.01)
    parser.add_argument("--vcf_region",type=str, default = "NA")
    parser.add_argument("--min_mapping_qual",type=int, default = 30)


    args = parser.parse_args()
    print args

    res_df = count_umis(args.input_vcf_file, args.bamfile, args)
    if res_df is not None:
        res_df.to_csv(args.output_file, sep = "\t", index = None)
    else:
        print "Warning: no variants detected"
    return res_df

# run count_umis_at_vcf_sites.py --input_vcf_file /data1/external_data/kees/Exac/ExAC.r0.3.1.sites.vep.vcf.gz --bamfile chr_split/pbmc33k_possorted_genome_bam.bam.chrom_1.rg.splitntrim.bam --output_file worst.txt --vcf_region 1:1-1000000

if __name__ == "__main__":
    res_df = main()

    print "Done!"
    

