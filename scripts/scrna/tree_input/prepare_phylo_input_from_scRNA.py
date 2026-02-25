import subprocess
import sys
import pandas as pd
from module.UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele
from module.read_file import handle_barcode
import os
from collections import defaultdict
import argparse
import pysam

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

def get_identifier_info_from_bam(identifier,bam_file,barcode_list,run_type="UMI"):
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
            CB=item.get_tag("CB").strip()
            UB=item.get_tag("UB").strip()
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
    
    used_barcode=site_barcode_UMI_dict.keys()
    
    allele_index_dict={"A":0,"T":1,"C":2,"G":3}
    count_info_dict=defaultdict(str)
    alt_info_dict=defaultdict(int)
    total_count=0
    alt_count=0
    if run_type=="UMI":
        for barcode_name in used_barcode:
            per_CB_list=[0,0,0,0]
            for UMI in site_barcode_UMI_dict[barcode_name]:
                count_dict=site_barcode_UMI_dict[barcode_name][UMI]["count"]
                quality_dict=site_barcode_UMI_dict[barcode_name][UMI]["quality"]
                phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
                candidate_allele,phred=get_most_candidate_allele(phred_dict,ref)
                per_CB_list[allele_index_dict[candidate_allele]]+=1
            
            count_info_dict[barcode_name]=str(per_CB_list[allele_index_dict[alt]])+"/"+str(sum(per_CB_list))
            total_count+=sum(per_CB_list)
            alt_info_dict[barcode_name]=per_CB_list[allele_index_dict[alt]]
            alt_count+=per_CB_list[allele_index_dict[alt]]
    
    elif run_type=="read":
        for barcode_name in used_barcode:
            per_read_list=[0,0,0,0]
            for UB in site_barcode_UMI_dict[CB].keys():
                for i,geno in zip([0,1,2,3],"ATCG"):
                    per_read_list[i]+=site_barcode_UMI_dict[CB][UB]["count"][geno]
            count_info_dict[barcode_name]=str(per_read_list[allele_index_dict[alt]])+"/"+str(sum(per_read_list))
            total_count+=sum(per_read_list)
            alt_info_dict[barcode_name]=per_read_list[allele_index_dict[alt]]
            alt_count+=per_read_list[allele_index_dict[alt]]
    
    count_values = [count_info_dict[key] if key in count_info_dict.keys() else '0/0' for key in barcode_list]
    posteria_values = [
        '1' if key in alt_info_dict.keys() and alt_info_dict[key] > 0 
        else '0' if key in alt_info_dict.keys() and alt_info_dict[key] == 0 
        else 'NA' 
        for key in barcode_list
    ]
    
    mut_likelihood_values = [
        '1' if key in alt_info_dict.keys() and alt_info_dict[key] > 0 
        else '0' if key in alt_info_dict.keys() and alt_info_dict[key] == 0 
        else 'NA' 
        for key in barcode_list
    ]
    
    nomut_likelihood_values = [
        '0' if key in alt_info_dict.keys() and alt_info_dict[key] > 0 
        else '1' if key in alt_info_dict.keys() and alt_info_dict[key] == 0 
        else 'NA' 
        for key in barcode_list
    ]
    
    
    # posteria_values = [ '1' if key in alt_info_dict.keys() and alt_info_dict[key] >0 else '0' for key in barcode_list]
    # mut_likelihood_values = [  '1' if key in alt_info_dict.keys() and alt_info_dict[key] >0 else '0' for key in barcode_list]
    # nomut_likelihood_values = [  '0' if key in alt_info_dict.keys() and alt_info_dict[key] >0 else '1' for key in barcode_list]
    # mut_likelihood_values = [ 'NA' for key in barcode_list]
    # nomut_likelihood_values = [ 'NA' for key in barcode_list]
    new_ref=str('"')+ref+str('"')
    
    return [chrom,pos,new_ref, alt, "1"], [alt_count, total_count], count_values, posteria_values, mut_likelihood_values, nomut_likelihood_values

def handle_posname(pos_name):
    sitem=pos_name.split("_")
    chrom=str(sitem[0])
    pos=int(sitem[1])
    ref=str(sitem[2])
    alt=str(sitem[3])
    
    return chrom, pos, ref, alt


def read_mutation_file(mutation_file):
    mutation_identifier_list=[]
    for line in open(mutation_file,"r"):
        sline=line.strip().split()
        if len(sline)==1:
            mutation_identifier_list.append(sline[0])
        elif len(sline)!=1:
            chrom=sline[0]; pos=sline[1]; ref= sline[2]; alt=sline[3]
            
            mutation_identifier_list.append("_".join([chrom,pos,ref,alt]))
    
    return mutation_identifier_list


# def get_ind_posteria_per_site(ind_genotype_file,mutation_identifier):
#     '''
#     mutation_identifier: chr13_94437452_G_A
#     ind_genotype_file:
#         #chrom  site    ID      germline        mutant  cluster spot_number     consensus_read_count    genotype        p_mosaic        Gi      vaf
#         GL000008.2      80303   .       C       G       bulk    2       0,0,1,1 mosaic  0.9999850002250079      0/1     0.5
#         GL000008.2      80321   .       T       A       bulk    2       1,1,0,0 mosaic  0.9999850002250079      0/1     0.5

#     '''
#     chrom,pos,ref,alt=handle_posname(mutation_identifier)
#     command="awk '$1==\"%s\" && $2==%s && $4==\"%s\" && $5==\"%s\" {print $8,$10}' %s " % (chrom, pos, ref, alt, ind_genotype_file)
 
#     result = subprocess.check_output(command,shell=True,text=True)
#     # print(result)
#     count,ind_posteria=result.strip().split(" ")[0],result.strip().split(" ")[1]
#     if ind_posteria=="" or count=="":
#         return [], []
    
#     count_list=count.split(",")
#     geno_index="ATCG".index(alt)
#     alt_count=count_list[geno_index]
#     all_count= sum(int(item) for item in count_list) 
#     new_ref=str('"')+ref+str('"')
#     total_count=str(alt_count)+"/"+str(all_count)

#     return [chrom,pos,new_ref, alt, ind_posteria], [total_count]


def tidy_func(sample,barcode_list, bam_file, run_type,mutation_identifier):
    ind_info, counts, count_values, posteria_values, mut_likelihood_values, nomut_likelihood_values =get_identifier_info_from_bam(mutation_identifier,bam_file,barcode_list,run_type)
    # out_list=[sample]+ ind_info + posteria_values + mut_likelihood_values + nomut_likelihood_values+ total_count+ count_values   
    # chrom,pos,new_ref, alt, ind_posteria , [posteria_values] , [mut_likelihood_values] , [nomut_likelihood_values], total_count, [count_values]
    return sample, ind_info , posteria_values , mut_likelihood_values , nomut_likelihood_values, counts, count_values
    # return out_list

def main():
    samples=str(args.samples).split(",")
    bam_files=str(args.bams).split(",")
    
    if args.barcode_files:
        barcode_files=str(args.barcode_files).split(",")
        barcode_lists=[]
        spot_number=0
        for i in barcode_files:
            barcode_dict=handle_barcode(i)
            spot_number+=len(list(barcode_dict.keys()))
            barcode_lists.append(list(barcode_dict.keys()))
    elif args.target_barcodes:
        barcode_lists=[]
        barcode_df=pd.read_csv(args.target_barcodes,sep="\t",header=None,names=["sample","barcode"])
        spot_number=barcode_df.shape[0]
        for sample in samples:
            barcode_lists.append(barcode_df[barcode_df["sample"] == sample]["barcode"].tolist())
    else:
        print("No barcode files or target barcodes provided.")
        sys.exit()
    
    
    # ind_genotype_file=args.ind_geno
    # sample=args.sample
    # bam_file=args.bam
    
    mutation_identifier_list=read_mutation_file(args.mutlist)
    out_list=[]
    
    out_name=args.outprefix+"_spot_c_"+str(spot_number)+".csv"
    out_file=open(out_name, "w")
    for mutation_identifier in mutation_identifier_list:
        posteria_values_list=[]
        mut_likelihood_values_list=[]
        nomut_likelihood_values_list=[]
        alt_counts=0; total_counts=0
        count_values_list=[]
        for sample,bam_file,barcode_list in zip(samples, bam_files,barcode_lists):
            _, ind_info , posteria_values , \
            mut_likelihood_values , nomut_likelihood_values, \
            counts, count_values = tidy_func(sample,barcode_list, bam_file,args.run_type,mutation_identifier)
            
            posteria_values_list=posteria_values_list+ posteria_values
            mut_likelihood_values_list=mut_likelihood_values_list+mut_likelihood_values
            nomut_likelihood_values_list=nomut_likelihood_values_list+nomut_likelihood_values
            alt_counts+=counts[0]; total_counts+=counts[1]
            count_values_list=count_values_list+count_values
        
        out_list.append(["combine"]+
                        ind_info+
                        posteria_values_list+
                        mut_likelihood_values_list+
                        nomut_likelihood_values_list+
                        [str(alt_counts)+"/"+str(total_counts)]+
                        count_values_list)
    
    for line in out_list:
        # print(line)
        write_info=",".join([str(k) for k in line if k != "" ])
        out_file.write(f'{write_info}\n')
    
    out_file.close()
    
    out_dir=os.path.dirname(args.outprefix)
    outsuffix=os.path.basename(args.outprefix)
    out_barcode_file=open(os.path.join(out_dir,outsuffix+"_scid_barcode.txt"),"w")
    out_barcode_file.write(f'scid_basedTree\n')
    for sample,barcode_list in zip(samples,barcode_lists):
        for barcode in barcode_list:
            out_barcode_file.write(f'{sample}_{barcode}\n')
    out_barcode_file.close()


## parameters
parser = argparse.ArgumentParser()

group = parser.add_mutually_exclusive_group(required=False)
group.add_argument("--barcode_files", help="barcode_file")
group.add_argument("--target_barcodes", help="one file contain the target barcodes, 1st col is sample, 2nd col is barcode")

parser.add_argument("--bams", required=True,help="bam file")
parser.add_argument("--mutlist", required=True,help="The file storing mutation list or from sf/phase result")
parser.add_argument("--outprefix", required=True,help="The output file prefix")
parser.add_argument("--samples", required=True,type=str,help="sample name")
parser.add_argument("--type", required=False,dest="run_type", default="UMI",choices=["UMI","read"],type=str,help="Do you want to use UMI or read")
args = parser.parse_args()


if __name__ == '__main__':
    main()

    
    