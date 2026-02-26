import argparse
from collections import defaultdict
from functools import reduce
from math import ceil, log10
import os
import subprocess
import pysam
import pandas as pd
from pybedtools import BedTool
import numpy as np
from UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele
from myutils import check_dir
from scipy.stats import binom
import statsmodels.stats.multitest as smm


def handle_cigar(ciagr_symbol):
    '''
    ## handel cigar
    # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
    # '76M1D33M139241N11M'
    # the 1st is symbol; and the 2nd is count
    # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
    output:
    set_cut may be a turple, contain the softclip information; pos_cut are the indel information: [(before_length, insert number)]
    '''
    seq_length_before = 0
    pos_length_before = 0

    seq_cut_start = None; seq_cut_end = None
    pos_cut=[]
    for cigars, i in zip(ciagr_symbol,range(1,len(ciagr_symbol)+1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol in [5,6,7,8]:
            # an api for handeling mapping issues "HP=X"
            print(ciagr_symbol)  ## LOG
        elif symbol in [0, 1, 4]:
            # measure the seq length 
            seq_length_before += count
            if symbol == 0:
                pos_length_before += count
            elif symbol == 4:
                # whether "S" is in this read
                if i == 1:
                    # whether the "S" is in the head or tail
                    seq_cut_start = seq_length_before
                elif i ==len(ciagr_symbol):
                    seq_cut_end = seq_length_before
                else:
                    print(ciagr_symbol) ## LOG
            elif symbol == 1:
                # whether the "I" is in the cigar
                pos_cut.append((pos_length_before,count))
        else:
            pass
    seq_cut = (seq_cut_start, seq_cut_end)
    return seq_cut, pos_cut


def handle_seq(seq, seq_cut):
    # only support one time for cut
    cut_seq=seq[seq_cut[0]:seq_cut[1]]
    return cut_seq


def handle_pos(pos_matrix,pos_cut):
    if len(pos_cut) == 0:
        cut_pos_matrix = pos_matrix
    elif len(pos_cut) == 1:
        times = pos_cut[0][1]
        pos = pos_cut[0][0]
        cut_pos_matrix = pos_matrix[0:pos] + [""] * times + pos_matrix[pos:]
    else:
        start = 0
        cut_pos_matrix = []
        for item in range(len(pos_cut)):
            pos = pos_cut[item][0]; times = pos_cut[item][1]
            cut_pos_matrix = cut_pos_matrix + pos_matrix[start:pos] + [""] * times
            start=pos
            #print(cut_pos_matrix)
        last_pos = pos_cut[-1][0]
        cut_pos_matrix = cut_pos_matrix + pos_matrix[last_pos:]
    return cut_pos_matrix


def handle_quality_matrix(mutation_in_cutseq_index,seq,cut_seq):
    if len(cut_seq[mutation_in_cutseq_index:]) >= len(cut_seq[:mutation_in_cutseq_index]):
        query_str = cut_seq[mutation_in_cutseq_index:]
        raw_index = seq.index(query_str)
    else:
        query_str = cut_seq[:mutation_in_cutseq_index]
        raw_index = seq.index(query_str) + len(query_str)
    return raw_index

def get_identifier_info_from_bam(identifier,bam_file,barcode_list,CBtag,UBtag,run_type="UMI"):
    # geno_index={"A":0,"T":1,"C":2,"G":3}
    bam_handle=pysam.AlignmentFile(bam_file,"r")
    chrom,pos,ref,alt=identifier.split("_")
    reads=bam_handle.fetch(chrom,int(pos)-1,int(pos))
    pos_index = int(pos)-1
    
    # all_reads_count=0
    # effective_DP=0
    # mapping_reads=0;skip_reads=0
    # barcode_low_quality=0   

    if barcode_list!=[]:
        filteration="yes"
    else:
        filteration="no"

    barcode_name=[]
    site_barcode_UMI_dict={}
    # counts=defaultdict(int)

    for item in reads:
        try:
            # a part of reads didn't have the information of "CB", bacause the "CR" didn't pass QC
            CB=item.get_tag(CBtag).strip()
            UB=item.get_tag(UBtag).strip()
        except:
                continue
  
        if filteration=="no" or (filteration=="yes" and CB in barcode_list):
            
            try:
                item.get_reference_positions().index(pos_index)
            except:
                # skip_reads += 1
                continue
            # mapping_reads += 1

            seq_cut, pos_cut = handle_cigar(item.cigar)
            cut_seq=handle_seq(item.seq, seq_cut)
            cut_pos=handle_pos(item.get_reference_positions(), pos_cut)

            if pos_index in cut_pos:
                # effective_DP += 1
                geno = cut_seq[cut_pos.index(pos_index)]
                if geno not in "ATCG":
                    continue
                
                raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
                quality=item.get_forward_qualities()[raw_index]

                UMI_name = str(UB)
                barcode_name = str(CB)

                if barcode_name not in site_barcode_UMI_dict.keys():
                    site_barcode_UMI_dict[barcode_name]=defaultdict(dict)

                if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                    site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                    site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}

                site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
                site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
    # print(site_barcode_UMI_dict)
    if run_type=="read":
        out_lists={}
        for CB in site_barcode_UMI_dict.keys():
            out_list=[0,0,0,0]
            for UB in site_barcode_UMI_dict[CB].keys():
                for i,geno in zip([0,1,2,3],"ATCG"):
                    out_list[i]+=site_barcode_UMI_dict[CB][UB]["count"][geno]
            out_lists[CB]=out_list
        return out_lists
    
    elif run_type=="check":
        out_lists={}
        for CB in site_barcode_UMI_dict.keys():
            for UB in site_barcode_UMI_dict[CB].keys():
                out_list=[0,0,0,0]
                for i,geno in zip([0,1,2,3],"ATCG"):
                    out_list[i]=site_barcode_UMI_dict[CB][UB]["count"][geno]
                out_lists[CB+":"+UB]=out_list
        return out_lists

    elif run_type=="UMI": 
        out_lists={}
        # consensus_read_count=defaultdict(int)
        for CB in site_barcode_UMI_dict.keys():
            out_list=[0,0,0,0]
            for UB in site_barcode_UMI_dict[CB].keys():
                count_dict=site_barcode_UMI_dict[CB][UB]["count"]
                quality_dict=site_barcode_UMI_dict[CB][UB]["quality"]
                phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
                candidate_allele,phred=get_most_candidate_allele(phred_dict,ref)
                out_list["ATCG".index(candidate_allele)]+=1
            out_lists[CB]=out_list

        return out_lists



def handle_barcode(barcode_file,pos=0,in_tissue_choose=0):
    '''
    input: (barcode,in_tissue_or_not,simplified_location_x,simplified_location_y,real_location_x,real_location_y) #in_tissue_or_not:0 means not in tissue, 1 means in tissue
    ACGCCTGACACGCGCT-1,0,0,0,323,308
    TACCGATCCAACACTT-1,0,1,1,334,326
    
    argrument:
    pos means which kind of position you want to use, if you given a int not equal to 0, the simplified location will be used.

    output:
    dict: {barcode: (in_tissue, pos1, pos2),...}
    '''
    barcode_dict={}
    f=open(barcode_file,"r")
    for line in f.readlines():
        s = line.strip().split(",")
        barcode=s[0]
        if barcode!="barcode":
            in_tissue=int(s[1])
            if pos==0:
                pos1= int(s[2]); pos2=int(s[3])
            else:
                pos1= int(s[4]); pos2=int(s[5])
            
            if in_tissue_choose==0 and in_tissue==1:
                barcode_dict[barcode]= (in_tissue, pos1, pos2)
            elif in_tissue_choose==1:
                barcode_dict[barcode]= (in_tissue, pos1, pos2)
            else:pass
        else: pass
    return barcode_dict

def main():
    file=open(args.mutation_list,"r")
    identifier_list=[i.strip() for i in file.readlines()]

    outdir=args.outdir
    check_dir(outdir)
    # bam_path="/storage/douyanmeiLab/yangzhirui/01.Data_download/01.Skin_Cancer/04.Analysis/01.cellranger_scrna_hg38"

    bed_file=os.path.join(outdir,"tmp.bed")
    run_shell="awk -F'_' '{print $1, $2, $2, $3, $4}' OFS='\t'  %s > %s" % (args.mutation_list, bed_file)
    result=subprocess.run(run_shell,shell=True)
    if result.returncode!=0:
        print(f'Something wrong when identifier tmp bed file for {bed_file}.')

    ase_file=args.ase_file
    # ase_file="/storage/douyanmeiLab/yangzhirui/01.Data_download/07.brain/04.Analysis/04.mutations/151673/candidate_ASE_sites.txt"
    ase_tmp_file=os.path.join(outdir,"ase.tmp.bed")
    run_shell = "sed '1d' %s |awk '{print $10, $11, $12, $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' - | sort -u | \
                 bedtools intersect -a - -b %s -wa > %s" % (ase_file,bed_file, ase_tmp_file)
    result=subprocess.run(run_shell,shell=True)
    if result.returncode!=0:
        print(f'Something wrong when ase tmp file for {ase_tmp_file}.')

    RNA_editing=args.RNA_editing
    RNA_editing_tmp_file=os.path.join(outdir,"RNA_editing.tmp.bed")
    run_shell = "bedtools intersect -a %s -b %s -wa > %s" % (RNA_editing,bed_file, RNA_editing_tmp_file)
    result=subprocess.run(run_shell,shell=True)
    if result.returncode!=0:
        print(f'Something wrong when RNA editing tmp file for {RNA_editing_tmp_file}.')
    

    if args.barcode_list!="" and os.path.exists(args.barcode_list):
            barcode_list=pd.read_csv(args.barcode_list,sep="\t",header=None,names=["barcode"])["barcode"].tolist()
    else:
        barcode_list=[]

    bam_file=args.bam
    run_type=args.run_type # read / check / UMI
    CBtag=args.CBtag
    UBtag=args.UBtag
    threshold=0.05
    out_file=open(outdir+"/"+"identifer_stat.txt","w")
    out_file.write(f'identifier\tase\tediting\n')

    for identifier in identifier_list:
        chrom,pos,ref,alt=identifier.split("_")

        out_info=get_identifier_info_from_bam(identifier,bam_file,barcode_list,CBtag,UBtag,run_type)
        countsum=0
        alt_count=0
        for onekey in out_info.keys():
            countsum += np.sum(out_info[onekey])
            alt_count += out_info[onekey]["ATCG".index(alt)] 
        # mosaic_vaf=int(alt_count)/countsum

        identifier_bed=BedTool("\t".join([chrom,pos,pos,ref,alt]),from_string=True)
        ase_tmp_bed=BedTool(ase_tmp_file)

        pos_intersect = identifier_bed.intersect(ase_tmp_bed, wb=True)

        if pos_intersect.count()==0:
            ase="Unknown"
        else:
            ase_list=[]
            directions=[]
            for line in pos_intersect:
                new_sline=line.fields
                print(new_sline)
                count1,count2=new_sline[12].split(",")
                germ_vaf=int(count1)/(int(count1)+int(count2))
                ase_p_value=float(new_sline[14])
                if germ_vaf<=0.35:
                    directions.append("ref")
                elif germ_vaf>=0.65:
                    directions.append("alt")
                else:
                    directions.append("balance")

                mosaic_p_value =  binom.cdf(int(alt_count), countsum, p=germ_vaf)
                ase_list.append(mosaic_p_value)

                # if ase_p_value>threshold:
                #     ase_list.append(0)
                # elif ase_p_value<=threshold and mosaic_p_value<=threshold:
                #     ase_list.append(1)
            
            ref_count = directions.count("ref")
            alt_count = directions.count("alt")
            balance_count = directions.count("balance")

            if ref_count > 0 and alt_count > 0:
                result = "unclear"  # 同时有 "ref" 和 "alt"
            elif ref_count > 0  and alt_count == 0:
                result = "ref"  # 只有 "ref" 和 "balance"
            elif alt_count > 0  and ref_count == 0:
                result = "alt"  # 只有 "alt" 和 "balance"
            elif balance_count == len(directions):
                result = "balance"  # 全部为 "balance"
            else:
                result = "unclear"  # 默认情况下为 unclear（如果符合其他情况）
            
            if result in ["unclear"]:
                ase="Unclear"
            elif result in ["balance"]:
                ase="False"
            else:
                _, p_adj, _, _ = smm.multipletests(ase_list, method='fdr_bh')
                
                if max(p_adj) >= threshold:
                    ase="True"
                else:
                    ase="False"

            # if result in ["ref","alt"] and sum(ase_list)>=1:
            #     ase="True"
            # elif result in ["balance"] or sum(ase_list)==0:
            #     ase="False"
            # else:
            #     ase="Unclear"

        RNA_editing_bed=BedTool(RNA_editing_tmp_file)
        editing_intersect = identifier_bed.intersect(RNA_editing_bed, wb=True)
        if editing_intersect.count()!=0:
            editing="True"
        else:
            editing="False"

        out_file.write(f'{identifier}\t{ase}\t{editing}\n')


## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--bam", required=True,help="bam_file")
parser.add_argument("--RNA_editing", required=True,help="RNA editing bed file")
parser.add_argument("--CBtag", required=False,default="CB",help="CB tag")
parser.add_argument("--UBtag", required=False,default="UB",help="UB tag")
parser.add_argument("--ase_file", required=True,help="ase file")
parser.add_argument("--barcode_list", required=False,default="",help="target barcodes file")
parser.add_argument("--run_type", required=False,choices=["UMI","read"],default="UMI",help="which type you want to run")
parser.add_argument("--mutation_list",required=False, help="identifier list: chr_pos_ref_alt")
parser.add_argument("--outdir","-o",required=False, default="IGV_Plot", help="output dir")

args = parser.parse_args()


if __name__ == '__main__':
    main()
