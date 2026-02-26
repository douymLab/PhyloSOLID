import pickle
import subprocess

import scipy
# from utils import  calculate_rbc_for_paired_wilcoxon, check_dir,combine_info_from_cigar, get_chr_size, get_indel_info, handle_p_value_log10,handle_posname,judge_pos_in_indel,do_wilicox_sum_test, round_to_nearest_bin, wilcoxon_with_rbc
import os
import pysam
from UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele, handle_cigar, handle_pos, handle_quality_matrix, handle_seq
from collections import Counter,defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc
import multiprocessing
from statistics import mean,median
from functools import partial
import argparse
import statsmodels.stats.multitest as smm

import sys
import os

# 添加你的项目根目录路径到 sys.path
sys.path.insert(0, "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/pre-classifier/scRNA/scripts/spatial-lineager")

from utils import calculate_rbc_for_paired_wilcoxon, check_dir, combine_info_from_cigar, get_chr_size, get_indel_info, handle_p_value_log10, handle_posname, judge_pos_in_indel, do_wilicox_sum_test, round_to_nearest_bin, wilcoxon_with_rbc, calculate_UMI_combine_phred, get_most_candidate_allele, handle_pos, handle_quality_matrix, handle_seq


def get_intergration_from_identifier_and_file(compare_pl_path,identifier,query_file,index_in_query="0 1 2 3"):
    """
    This function is used to get the intergration result from identifier and one special query file.
    Note: The top 4 columns of query_file must be "chrom pos ref alt"
    
    """
    perl_script=compare_pl_path
    identifier_line="\t".join(identifier.strip().split("_"))

    command = f'echo -e "{identifier_line}" | perl {perl_script} - {query_file} {index_in_query} | sort -u'
    try:
        result=subprocess.check_output(command,text=True,shell=True)
        return result
    except:
        print(f"Something wrong when run the command: {command}")
        return ""


def refine_mean(in_list):
    try:
        return mean(in_list)
    except:
        print("wrong in refine mean for", in_list)
        return "NA"


def refine_median(in_list):
    try:
        return median(in_list)
    except:
        print("wrong in refine median for", in_list)
        return "NA"
        

def check_UMIconsistence_for_each_geno(count_dict,threshold=1):
    '''
    This function is used to count the consistence or not for each geno and each dict
    Version1: we want to contaion those info: A:9,T:1. Both geno A and T will be counted as 1 UMI inconsistence
    '''
    UMI_DP=sum(count_dict.values())
    if UMI_DP>=threshold:
        norm_count=[count_dict[geno]/UMI_DP for geno in "ATCG"]
        return norm_count
    else:
        return []
    


def robust_standardization(data):
    """使用中位数和IQR对数据进行标准化"""
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    standardized_data = (data - median) / iqr
    return standardized_data

def permutation_test(list1, list2, num_permutations=1000):
    """Permutation Test 用于评估两组数据的差异"""
    # 计算原始的均值差异
    observed_diff = np.mean(list1) - np.mean(list2)
    
    # 合并两组数据进行置换
    combined_data = np.concatenate([list1, list2])
    count = 0
    
    for _ in range(num_permutations):
        # 打乱数据并重新分成两组
        np.random.shuffle(combined_data)
        new_list1 = combined_data[:len(list1)]
        new_list2 = combined_data[len(list1):]
        
        # 计算新的均值差异
        permuted_diff = np.mean(new_list1) - np.mean(new_list2)
        
        # 计算置换差异大于等于原始差异的次数
        if np.abs(permuted_diff) >= np.abs(observed_diff):
            count += 1
    
    # 计算 p 值
    p_value = count / num_permutations
    return observed_diff, p_value


def handle_bam_file(bam_file,chrom,pos,ref,alt,run_type,outdir,bins,readLen=120):
    exist_CB=[]
    if "," in ref:
        one_ref=ref[0]
    else:
        one_ref=ref

    result_dict={"A":defaultdict(list), "T":defaultdict(list), "C":defaultdict(list), "G":defaultdict(list), "del":defaultdict(list)}
    for geno in "ATCG":
        result_dict[geno]["dp"]=0
        result_dict[geno]["dp_consensus"]=0
        result_dict[geno]["reverse_dp"]=0
        result_dict[geno]["forward_dp"]=0
        result_dict[geno]["edist"]=[0]*readLen
        result_dict[geno]["GenoSpotNum"]=0

    dp=0
    in_bam_read=pysam.AlignmentFile(bam_file, "rb") # , reference_filename=ref_fasta)
    pos_index = pos-1
    barcode_name=[]
    site_barcode_UMI_dict={}
    for read in in_bam_read.fetch(chrom, pos-1, pos):
        if run_type=="visium":
            try:
                CB=read.get_tag("CB").strip()
                UB=read.get_tag("UB").strip()

                barcode_name=str(CB)
                UMI_name=str(UB)
            except:
                continue
        elif run_type=="stereo":
            try:
                Cx_raw=int(read.get_tag("Cx"))
                Cy_raw=int(read.get_tag("Cy"))
                if bins !=1:
                    Cx=round_to_nearest_bin(Cx_raw,bins)
                    Cy=round_to_nearest_bin(Cy_raw,bins)
                else:
                    Cx=Cx_raw
                    Cy=Cy_raw
                UR=read.get_tag("UR").strip()

                barcode_name=str(Cx)+"_"+str(Cy)
                UMI_name=str(Cx_raw)+"_"+str(Cy_raw)+"_"+str(UR)
            except:
                continue

        elif run_type=="ST":
            try:
                CB=str(read.get_tag("B0"))
                UB=str(read.get_tag("B3"))
                barcode_name=str(CB)
                UMI_name=str(UB)
            except:
                continue
        else:
            # print("type",run_type)
            continue

        pos_index=pos-1
        seq_soft_cut, ins_info, del_info, seq_hard_clip = combine_info_from_cigar(read.cigar)
        cut_seq=handle_seq(read.seq, seq_soft_cut)
        cut_pos=handle_pos(read.get_reference_positions(), ins_info)
        indel_pos_list=judge_pos_in_indel(ins_info,del_info,read.get_reference_positions())

        if pos_index in cut_pos or pos_index in indel_pos_list:
            dp+=1
            if pos_index in indel_pos_list:
                geno="del"
                result_dict[geno]["is_indel"].append(1)
                result_dict[geno]["baseq"]=[]

            else:
                edist=cut_pos.index(pos_index)  
                geno = cut_seq[edist]  
                if geno not in "ATCG":
                    continue

                epos=edist/len(cut_pos)
                result_dict[geno]["epos"].append(epos)
                result_dict[geno]["edist"][edist]+=1

                result_dict["is_indel"]=[]
                raw_index = handle_quality_matrix(cut_pos.index(pos_index),read.seq,cut_seq)
                quality=read.get_forward_qualities()[raw_index]
                result_dict[geno]["baseq"].append(quality)
                # if result_dict[geno]["dp"]!=[]:
                result_dict[geno]["dp"]+=1
                # else:
                #     result_dict[geno]["dp"]=0
            # effective_DP += 1
        
            #number_mismatch; is_reverse; mapping_quality
                number_mismatch=read.get_tag("nM"); result_dict[geno]["number_mismatch"].append(number_mismatch)
                is_reverse=read.is_reverse; result_dict[geno]["is_reverse"].append(is_reverse)
                map_q=read.mapq; result_dict[geno]["map_q"].append(map_q)
                
                number_mapper=read.get_tag("NH"); result_dict[geno]["number_mapper"].append(number_mapper)

                #soft_clip_length and hard_clip_length
                left_softclip=0 if seq_soft_cut[0]==None else seq_soft_cut[0]
                right_softclip=0 if seq_soft_cut[1]==None else len(read.seq)-seq_soft_cut[1]
                softclip_length=left_softclip+right_softclip
                result_dict[geno]["left_softclip"].append(left_softclip)
                result_dict[geno]["right_softclip"].append(right_softclip)
                result_dict[geno]["softclip_length"].append(softclip_length)

                left_hardclip,right_hardclip=seq_hard_clip[0],seq_hard_clip[1]
                hardclip_length=left_hardclip+right_hardclip
                result_dict[geno]["left_hardclip"].append(left_hardclip)
                result_dict[geno]["right_hardclip"].append(right_hardclip)
                result_dict[geno]["hardclip_length"].append(hardclip_length)

                #indel information, indel number, indel length, indel distance
                ins_num,ins_length,ins_distance=get_indel_info(ins_info,read.get_reference_positions().index(pos_index))
                del_num,del_length,del_distance=get_indel_info(del_info,read.get_reference_positions().index(pos_index))
                result_dict[geno]["ind_num"].append(ins_num+del_num)
                result_dict[geno]["ins_num"].append(ins_num)
                if ins_num==0:
                    result_dict[geno]["ins_length"].append(["no"]); result_dict[geno]["ins_distance"].append(["no"])
                else: #the ins_length and ins_distance are list format
                    result_dict[geno]["ins_length"].append(ins_length) ## append a list
                    result_dict[geno]["ins_distance"].append(ins_distance) ## append a list
                
                result_dict[geno]["del_num"].append(del_num)
                if del_num==0:
                    result_dict[geno]["del_length"].append(["no"]); result_dict[geno]["del_distance"].append(["no"])
                else: # the del_length and del_distance are list format
                    result_dict[geno]["del_length"].append(del_length)
                    result_dict[geno]["del_distance"].append(del_distance)

                # querypos(querypos_p): the distance between pos and read start (doubt: the more far away from 1st seq pos, the lower quality may have), 
                # seqpos_p cycling length, related with strand (note: next_reference_start is only work for PE); 
                # for visium, all reads are read2, so seqpos may same as the len(querypos)
                # left pos: mapping position for the reference start; 
                left_boundary=edist+left_softclip+left_hardclip
                right_boundary=len(cut_pos)-edist + right_softclip + right_hardclip
                result_dict[geno]["left_read_edist"].append(edist)
                result_dict[geno]["right_read_edist"].append(len(cut_pos)-edist)

                left_boundary_remove_clip=edist
                right_boundary_remove_clip=len(cut_pos)-edist
                result_dict[geno]["querypos"].append(left_boundary)
                result_dict[geno]["seqpos"].append(right_boundary)
                if is_reverse in [True,"TRUE","true","True"]:
                    distance_to_end=right_boundary/readLen
                    result_dict[geno]["reverse_dp"]+=1
                    distance_to_end_remove_clip=right_boundary_remove_clip/len(cut_pos)
                    distance_to_end_remove_clip_save=right_boundary_remove_clip
                else:
                    distance_to_end=left_boundary/readLen
                    result_dict[geno]["forward_dp"]+=1
                    distance_to_end_remove_clip=left_boundary_remove_clip/len(cut_pos)
                    distance_to_end_remove_clip_save=left_boundary_remove_clip
                # print(geno,distance_to_end_remove_clip)
                result_dict[geno]["distance_to_end"].append(distance_to_end)
                result_dict[geno]["distance_to_end_remove_clip"].append(distance_to_end_remove_clip)

                leftpos_p=read.reference_start
                rightpos_p=read.reference_end # same as leftpo, can be deleted 
                result_dict[geno]["leftpos_p"].append(leftpos_p)
                result_dict[geno]["rightpos_p"].append(rightpos_p)

                #baseq1b
                if pos_index+1 in cut_pos:
                    baseq1b=read.get_forward_qualities()[raw_index+1]
                else:
                    baseq1b=""
                result_dict[geno]["baseq1b"].append(baseq1b)
                # print(read)
                #gene information
                try:
                    result_dict[geno]["GeneID_list"].append(read.get_tag("GX"))
                except:
                    result_dict[geno]["GeneID_list"].append("no")
                try:
                    result_dict[geno]["GeneName_list"].append(read.get_tag("GN"))
                except:
                    result_dict[geno]["GeneName_list"].append("no")
                try:
                    #'ENST00000301072,+1576,120M;ENST00000541364,+1539,120M;ENST00000552448,+1650,120M;ENST00000639419,+923,120M')
                    for item in read.get_tag("TX").split(";"):
                        transcript_id,_,_=item.split(",")
                        result_dict[geno]["TransID_list"].append(transcript_id)
                except:
                    result_dict[geno]["TransID_list"].append("no")

                if barcode_name not in site_barcode_UMI_dict.keys():
                    site_barcode_UMI_dict[barcode_name]=defaultdict(dict)

                if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                    site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                    site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
                    # site_barcode_UMI_dict[barcode_name][UMI_name]["context"]=[]
                    site_barcode_UMI_dict[barcode_name][UMI_name]["end"]=[]
                    site_barcode_UMI_dict[barcode_name][UMI_name]["end_remove_clip"]=[]
                    

                site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
                site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
                # site_barcode_UMI_dict[barcode_name][UMI_name]["context"].append(cut_seq[max(0,read_index-4):min(read_index+5,len(cut_seq))])
                site_barcode_UMI_dict[barcode_name][UMI_name]["end"].append(distance_to_end)
                site_barcode_UMI_dict[barcode_name][UMI_name]["end_remove_clip"].append(distance_to_end_remove_clip)
                
    for barcode in site_barcode_UMI_dict.keys():
        read_have_alt=False
        read_number_per_spot=0
        UMI_number_per_spot=0
        alt_UMI_number_per_spot=0
        UMI_count_by_allele=[0,0,0,0]
        # UMI_dp+=len(site_barcode_UMI_dict[barcode].keys())
        for UMI in site_barcode_UMI_dict[barcode]:             
            count_dict=site_barcode_UMI_dict[barcode][UMI]["count"]
            quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
            phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
            candidate_allele,phred=get_most_candidate_allele(phred_dict,one_ref)
            result_dict[candidate_allele]["dp_consensus"]+=1
            UMI_count_by_allele["ATCG".index(candidate_allele)]+=1
            threshold=1
            norm_count=check_UMIconsistence_for_each_geno(count_dict,threshold)
            norm_count_remove_single_read=check_UMIconsistence_for_each_geno(count_dict,2)

            if norm_count!=[]:
                for geno,prop in zip("ATCG",norm_count):
                    result_dict[geno]["UMI_consistence_prop"].append(prop)

            if norm_count_remove_single_read!=[]:
                for geno,prop in zip("ATCG",norm_count_remove_single_read):
                    result_dict[geno]["UMI_consistence_prop_remove_single_read"].append(prop)

            for geno in "ATCG":
                if site_barcode_UMI_dict[barcode][UMI]["count"][geno]!=0:
                    # print([count_dict["A"],count_dict["T"],count_dict["C"],count_dict["G"]])
                    result_dict[geno]["read_number_per_UMI"].append(count_dict[geno])
                    read_number_per_spot+=site_barcode_UMI_dict[barcode][UMI]["count"][geno]
                    
            UMI_number_per_spot+=1
            end=np.median(site_barcode_UMI_dict[barcode][UMI]["end"])
            end_remove_clip=np.median(site_barcode_UMI_dict[barcode][UMI]["end_remove_clip"])

            if candidate_allele==alt:
                read_have_alt=True
                alt_UMI_number_per_spot+=1
        
        if read_have_alt==True:
            # print(barcode)
            result_dict[alt]["GenoSpotNum"]+=1 
            result_dict[alt]["total_read_number_per_spot"].append(read_number_per_spot)
            result_dict[alt]["total_UMI_number_per_spot"].append(UMI_number_per_spot)
            result_dict[alt]["UMI_end"].append(end)
            result_dict[alt]["UMI_end_remove_clip"].append(end_remove_clip)
            result_dict[alt]["vaf_spot"].append(alt_UMI_number_per_spot/UMI_number_per_spot)

        else:
            result_dict[one_ref]["GenoSpotNum"]+=1
            result_dict[one_ref]["total_read_number_per_spot"].append(read_number_per_spot)
            result_dict[one_ref]["total_UMI_number_per_spot"].append(UMI_number_per_spot)
            result_dict[one_ref]["UMI_end"].append(end)
            result_dict[one_ref]["UMI_end_remove_clip"].append(end_remove_clip)
        
        for geno,count in zip("ATCG",UMI_count_by_allele):
            result_dict[geno]["UMI_number_per_spot"].append(count)

    def simple_to_get_list(read_info_dict,ref_allele,alt_allele,var):
        # print(ref_allele,alt_allele, ref_allele.split(","))
        return_ref_list=[float(k) for allele in ref_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
        return_alt_list=[float(k) for allele in alt_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
        return return_ref_list,return_alt_list

    get_list=partial(simple_to_get_list,result_dict,ref,alt)

    odds_eps=0
    ### new features
    ref_mapq,alt_mapq=get_list("map_q")
    ref_mapq_mean=refine_mean(ref_mapq)
    alt_mapq_mean=refine_mean(ref_mapq)
    mapq_mean=refine_mean(ref_mapq+alt_mapq)
    mapq_s,mapq_p,mapq_rbc=wilcoxon_with_rbc(ref_mapq,alt_mapq,alternative="greater")

    ref_baseq,alt_baseq=get_list("baseq")
    baseq_s,baseq_p,baseq_rbc=wilcoxon_with_rbc(ref_baseq,alt_baseq,alternative="greater")

    ref_baseq1b,alt_baseq1b=get_list("baseq1b")
    # ref_baseq1b_s,ref_baseq1b_p,ref_baseq1b_rbc=wilcoxon_with_rbc(ref_baseq,ref_baseq1b,alternative="greater")
    # alt_baseq1b_s,alt_baseq1b_p,alt_baseq1b_rbc=wilcoxon_with_rbc(alt_baseq,alt_baseq1b,alternative="greater")
    ref_baseq1b_s,ref_baseq1b_p,ref_baseq1b_rbc=wilcoxon_with_rbc(ref_baseq1b,ref_baseq,alternative="greater")
    alt_baseq1b_s,alt_baseq1b_p,alt_baseq1b_rbc=wilcoxon_with_rbc(alt_baseq1b,alt_baseq,alternative="greater")

    ref_querypos,alt_querypos=get_list("querypos")
    querypos_s,querypos_p,querypos_rbc=wilcoxon_with_rbc(ref_querypos,alt_querypos)

    ref_leftpos,alt_leftpos=get_list("leftpos_p")
    leftpos_s,leftpos_p,leftpos_rbc=wilcoxon_with_rbc(ref_leftpos,alt_leftpos)

    ref_seqpos,alt_seqpos=get_list("seqpos")
    seqpos_s,seqpos_p,seqpos_rbc=wilcoxon_with_rbc(ref_seqpos,alt_seqpos)

    ref_distance_to_end,alt_distance_to_end=get_list("distance_to_end")
    distance_to_end_s,distance_to_end_p,distance_to_end_rbc=wilcoxon_with_rbc(ref_distance_to_end,alt_distance_to_end,alternative="less")

    ref_distance_to_end_by_UMI,alt_distance_to_end_by_UMI=get_list("UMI_end")
    UMI_end_s,UMI_end_p,UMI_end_rbc=wilcoxon_with_rbc(ref_distance_to_end_by_UMI,alt_distance_to_end_by_UMI,alternative="less")

    ref_distance_to_end_remove_clip_by_UMI,alt_distance_to_end_remove_clip_by_UMI=get_list("UMI_end_remove_clip")
    UMI_end_remove_clip_s,UMI_end_remove_clip_p,UMI_end_remove_clip_rbc=wilcoxon_with_rbc(ref_distance_to_end_remove_clip_by_UMI,alt_distance_to_end_remove_clip_by_UMI,alternative="less")

    ref_mismatches, alt_mismatches=get_list("number_mismatch")
    refine_alt_mismatches=[int(i)-1 if i!=0 else 0 for i in alt_mismatches]
    mismatches_s,mismatches_p,mismatches_rbc=wilcoxon_with_rbc(ref_mismatches,refine_alt_mismatches)

    ref_softclip_length,alt_softclip_length=get_list("softclip_length")
    softclip_length_s,softclip_length_p,softclip_length_rbc=wilcoxon_with_rbc(ref_softclip_length,alt_softclip_length)

    all_softclip=ref_softclip_length+alt_softclip_length
    ref_softclip_prop=len([x for x in ref_softclip_length if x > 10])/len(ref_softclip_length) if len(ref_softclip_length)!=0 else "NA"
    alt_softclip_prop=len([x for x in alt_softclip_length if x > 10])/len(alt_softclip_length) if len(alt_softclip_length)!=0 else "NA"
    softclip_prop=len([x for x in all_softclip if x > 10])/len(all_softclip) if len(all_softclip)!=0 else "NA"
    # self.ref_softclip_length,self.alt_softclip_length=ref_softclip_length,alt_softclip_length
    ref_soft_count=len([x for x in ref_softclip_length if x > 10]); ref_no_soft_count=len(ref_softclip_length)-ref_soft_count
    alt_soft_count=len([x for x in alt_softclip_length if x > 10]); alt_no_soft_count=len(alt_softclip_length)-alt_soft_count
    softclip_prop_odds, softclip_prop_p=scipy.stats.fisher_exact([[alt_soft_count, alt_no_soft_count ],[ref_soft_count, ref_no_soft_count]])
        
    consensus_ref_allele_count=sum([float(result_dict[allele]["dp_consensus"]) for allele in ref.split(",") ])
    consensus_alt_allele_count=sum([float(result_dict[allele]["dp_consensus"]) for allele in alt.split(",") ])
    other_count_dict={base:int(result_dict[base]["dp_consensus"]) for base in "ATCG" if base not in ref+alt}
    consensus_alt2_allele_count=max(list(other_count_dict.values()),default=0)

    ref_allele_count=sum([float(result_dict[allele]["dp"]) for allele in ref.split(",") ])
    alt_allele_count=sum([float(result_dict[allele]["dp"]) for allele in alt.split(",") ])
    other_count_dict={base:int(result_dict[base]["dp"]) for base in "ATCG" if base not in ref+alt}
    alt2_allele_count=max(list(other_count_dict.values()),default=0)
    
    ref_read_number_perUMI,alt_read_number_perUMI=get_list("read_number_per_UMI")
    read_number_s,read_number_p,read_number_rbc=wilcoxon_with_rbc(ref_read_number_perUMI,alt_read_number_perUMI)

    ref_ind_num,alt_ind_num=get_list("ind_num")
    ref_read_with_indel = len([i for i in ref_ind_num if i > 0]); ref_read_without_indel = len(ref_ind_num)-ref_read_with_indel
    alt_read_with_indel = len([i for i in alt_ind_num if i > 0]); alt_read_without_indel = len(alt_ind_num)-alt_read_with_indel
    reads_with_indel_odds, reads_with_indel_p = scipy.stats.fisher_exact([[alt_read_with_indel+odds_eps, alt_read_without_indel+odds_eps ],[ref_read_with_indel+odds_eps, ref_read_without_indel+odds_eps]])

    ref_mappers, alt_mappers=get_list("number_mapper")
    refine_ref_mappers_uniq = len([i for i in ref_mappers if i == 1]); refine_ref_mappers_multi=len(ref_mappers)-refine_ref_mappers_uniq
    refine_alt_mappers_uniq = len([i for i in alt_mappers if i == 1]); refine_alt_mappers_multi=len(alt_mappers)-refine_alt_mappers_uniq
    multi_mapper_odds, multi_mapper_p = scipy.stats.fisher_exact([[refine_alt_mappers_multi+odds_eps, refine_alt_mappers_uniq+odds_eps ],[refine_ref_mappers_multi+odds_eps, refine_ref_mappers_uniq+odds_eps]])
    all_mappers=len(ref_mappers)+len(alt_mappers)
    multi_map_prop=(refine_ref_mappers_multi+refine_alt_mappers_multi)/all_mappers if all_mappers!=0 else "NA"
    # ref_multi_map_prop=refine_ref_mappers_multi/len(ref_mappers) if len(ref_mappers)!=0 else "NA"
    # alt_multi_map_prop=refine_alt_mappers_multi/len(alt_mappers) if len(alt_mappers)!=0 else "NA"
    
    ref_is_reverse, alt_is_reverse=get_list("is_reverse")
    refine_ref_is_reverse = len([i for i in ref_is_reverse if i == 1]); refine_ref_is_forward = len(ref_is_reverse)-refine_ref_is_reverse
    refine_alt_is_reverse = len([i for i in alt_is_reverse if i == 1]); refine_alt_is_forward = len(alt_is_reverse)-refine_alt_is_reverse
    strand_bias_odds, strand_bias_p = scipy.stats.fisher_exact([[refine_alt_is_reverse+odds_eps, refine_alt_is_forward+odds_eps ],[refine_ref_is_reverse+odds_eps, refine_ref_is_forward+odds_eps]])

    alt_UMI_consistence_remove_single_read=0; total_UMI_contain_alt_remove_single_read=0
    alt_consistence_hard_remove_single_read=0; alt_consistence_soft_remove_single_read=0
    for allele in alt.split(","):
        for k in result_dict[allele]["UMI_consistence_prop_remove_single_read"]:
            if k !=0.0:
                alt_UMI_consistence_remove_single_read+=k
                total_UMI_contain_alt_remove_single_read+=1
                if k==1.0:
                    alt_consistence_hard_remove_single_read+=1
                    alt_consistence_soft_remove_single_read+=1
                elif k>=0.75:
                    alt_consistence_soft_remove_single_read+=1
        try:
            alt_UMI_avg_consistence_remove_single_read=alt_UMI_consistence_remove_single_read/total_UMI_contain_alt_remove_single_read
            alt_consistent_UMI_prop_relaxed_remove_single_read=alt_consistence_soft_remove_single_read/total_UMI_contain_alt_remove_single_read
            alt_consistent_UMI_prop_strict_remove_single_read=alt_consistence_hard_remove_single_read/total_UMI_contain_alt_remove_single_read
    
        except:
            alt_UMI_avg_consistence_remove_single_read="no"
            alt_consistent_UMI_prop_relaxed_remove_single_read="no"
            alt_consistent_UMI_prop_strict_remove_single_read="no"
            pass
    
    alt_umi_number=[float(k) for allele in alt.split(",") for k in result_dict[allele]["UMI_number_per_spot"]]
    # print(alt_umi_number)
    umi_depth=[a+t+c+g for a,t,c,g in zip(result_dict["A"]["UMI_number_per_spot"], \
                                          result_dict["T"]["UMI_number_per_spot"], \
                                          result_dict["C"]["UMI_number_per_spot"], \
                                          result_dict["G"]["UMI_number_per_spot"])]
    alt_vs_total_dp_paired_wilcoxon_rbc=calculate_rbc_for_paired_wilcoxon(umi_depth,alt_umi_number)


    # ref_spot_dp_list, alt_spot_dp_list=get_list("total_UMI_number_per_spot")
    # UMI_number_per_spot_s,UMI_number_per_spot_p = wilcoxon_with_rbc(ref_spot_dp_list, alt_spot_dp_list,method="two-sided",type="list")
    # ref_UMI_number_per_spot_median,alt_UMI_number_per_spot_median = np.median(ref_spot_dp_list),np.median(alt_spot_dp_list)
    # ref_UMI_number_per_spot_max,alt_UMI_number_per_spot_max = max([0]+ref_spot_dp_list),max([0]+alt_spot_dp_list)

    vaf_spot_mean=refine_mean(result_dict[alt]["vaf_spot"])
    vaf_spot_median=refine_median(result_dict[alt]["vaf_spot"])
    vaf_values = np.array(result_dict[alt]["vaf_spot"])  
    only_alt_mutant_spot_prop = (vaf_values == 1).mean()

    p_list=[baseq_p,ref_baseq1b_p,alt_baseq1b_p,querypos_p,leftpos_p,seqpos_p,distance_to_end_p,UMI_end_p,UMI_end_remove_clip_p,mismatches_p,mapq_p,read_number_p,softclip_length_p]
    _, p_adj_list, _, _ = smm.multipletests(p_list, method='fdr_bh')
    baseq_p_adj,ref_baseq1b_p_adj,alt_baseq1b_p_adj,querypos_p_adj,leftpos_p_adj,seqpos_p_adj,distance_to_end_p_adj,UMI_end_p_adj,UMI_end_remove_clip_p_adj,mismatches_p_adj,mapq_p_adj,read_number_p_adj,softclip_length_p_adj=p_adj_list

    return baseq_s,baseq_p,baseq_rbc, ref_baseq1b_s,ref_baseq1b_p,ref_baseq1b_rbc, alt_baseq1b_s,alt_baseq1b_p,alt_baseq1b_rbc, \
                querypos_s,querypos_p,querypos_rbc,leftpos_s,leftpos_p,leftpos_rbc, seqpos_s, seqpos_p, seqpos_rbc, \
                distance_to_end_s,distance_to_end_p,distance_to_end_rbc, UMI_end_s,UMI_end_p,UMI_end_rbc, \
                UMI_end_remove_clip_s,UMI_end_remove_clip_p,UMI_end_remove_clip_rbc, \
                mismatches_s,mismatches_p,mismatches_rbc,mapq_s,mapq_p,mapq_rbc, mapq_mean, \
                alt_UMI_avg_consistence_remove_single_read,alt_consistent_UMI_prop_relaxed_remove_single_read,alt_consistent_UMI_prop_strict_remove_single_read, \
                read_number_s,read_number_p, read_number_rbc, strand_bias_odds, strand_bias_p,multi_mapper_odds, multi_mapper_p, \
                softclip_length_s,softclip_length_p,softclip_length_rbc, \
                ref_allele_count, consensus_ref_allele_count, alt_allele_count, consensus_alt_allele_count, \
                alt2_allele_count,consensus_alt2_allele_count,alt_vs_total_dp_paired_wilcoxon_rbc,reads_with_indel_odds, reads_with_indel_p, \
                ref_softclip_prop,alt_softclip_prop,softclip_prop,softclip_prop_odds,softclip_prop_p,ref_mapq_mean,alt_mapq_mean,multi_map_prop, \
                vaf_spot_mean,vaf_spot_median,only_alt_mutant_spot_prop, \
                baseq_p_adj,ref_baseq1b_p_adj,alt_baseq1b_p_adj,querypos_p_adj,leftpos_p_adj,seqpos_p_adj,distance_to_end_p_adj,UMI_end_p_adj,UMI_end_remove_clip_p_adj,mismatches_p_adj,mapq_p_adj,read_number_p_adj,softclip_length_p_adj, \
                {}


def handel_identifier(bam_file,run_type,readLen,outdir,bins,prior,identifier):
    # print("running",identifier)
    current_directory = os.path.dirname(os.path.abspath(__file__))
    compare_pl_path=os.path.join(os.path.dirname(current_directory),"others/compare_files.pl")
    chrom,pos,ref,alt=identifier.split("_")
    only_pos_identifier="\t".join([chrom,str(pos),chrom,str(pos)])

    features = _global_features[identifier]

    # columns = ['identifier','s1','p1','s2','p2','s3','p3','chi2_stat_read',
    #             'chi2_p_read','chi2_stat_UMI','chi2_p_UMI','mean_distance_to_end',
    #             'median_distance_to_end','mean_distance_to_end_remove_clip',
    #             'median_distance_to_end_remove_clip']
    baseq_s,baseq_p,baseq_rbc, ref_baseq1b_s,ref_baseq1b_p,ref_baseq1b_rbc, alt_baseq1b_s,alt_baseq1b_p,alt_baseq1b_rbc, \
        querypos_s,querypos_p,querypos_rbc,leftpos_s,leftpos_p,leftpos_rbc, seqpos_s, seqpos_p, seqpos_rbc, \
        distance_to_end_s,distance_to_end_p,distance_to_end_rbc, UMI_end_s,UMI_end_p,UMI_end_rbc, \
        UMI_end_remove_clip_s,UMI_end_remove_clip_p,UMI_end_remove_clip_rbc, \
        mismatches_s,mismatches_p,mismatches_rbc,mapq_s,mapq_p,mapq_rbc, mapq_mean, \
        alt_UMI_avg_consistence_remove_single_read,alt_consistent_UMI_prop_relaxed_remove_single_read,alt_consistent_UMI_prop_strict_remove_single_read, \
        read_number_s,read_number_p, read_number_rbc, strand_bias_odds, strand_bias_p,multi_mapper_odds, multi_mapper_p, \
        softclip_length_s,softclip_length_p,softclip_length_rbc, \
        ref_allele_count, consensus_ref_allele_count, alt_allele_count, consensus_alt_allele_count, \
        alt2_allele_count,consensus_alt2_allele_count,alt_vs_total_dp_paired_wilcoxon_rbc,reads_with_indel_odds, reads_with_indel_p, \
        ref_softclip_prop,alt_softclip_prop,softclip_prop,softclip_prop_odds,softclip_prop_p,ref_mapq_mean,alt_mapq_mean,multi_map_prop, \
        vaf_spot_mean,vaf_spot_median,only_alt_mutant_spot_prop, \
        baseq_p_adj,ref_baseq1b_p_adj,alt_baseq1b_p_adj,querypos_p_adj,leftpos_p_adj,seqpos_p_adj,distance_to_end_p_adj,UMI_end_p_adj,UMI_end_remove_clip_p_adj,mismatches_p_adj,mapq_p_adj,read_number_p_adj,softclip_length_p_adj, \
        result_dict =handle_bam_file(bam_file,chrom,int(pos),ref,alt,run_type,outdir,bins,readLen)
        
    base_data = {
        'baseq_s': baseq_s,
        'baseq_p': handle_p_value_log10(baseq_p),
        'baseq_rbc': baseq_rbc,

        'ref_baseq1b_s': ref_baseq1b_s,
        'ref_baseq1b_p': handle_p_value_log10(ref_baseq1b_p),
        'ref_baseq1b_rbc': ref_baseq1b_rbc,

        'alt_baseq1b_s': alt_baseq1b_s,
        'alt_baseq1b_p': handle_p_value_log10(alt_baseq1b_p),
        'alt_baseq1b_rbc': alt_baseq1b_rbc,

        'querypos_s': querypos_s,
        'querypos_p': handle_p_value_log10(querypos_p),
        'querypos_rbc': querypos_rbc,

        'leftpos_s': leftpos_s,
        'leftpos_p': handle_p_value_log10(leftpos_p),
        'leftpos_rbc': leftpos_rbc,

        'seqpos_s': seqpos_s,
        'seqpos_p': handle_p_value_log10(seqpos_p),
        'seqpos_rbc': seqpos_rbc,

        'distance_to_end_s': distance_to_end_s,
        'distance_to_end_p': handle_p_value_log10(distance_to_end_p),
        'distance_to_end_rbc': distance_to_end_rbc,

        'UMI_end_s': UMI_end_s,
        'UMI_end_p': handle_p_value_log10(UMI_end_p),
        'UMI_end_rbc': UMI_end_rbc,

        'mismatches_s': mismatches_s,
        'mismatches_p': handle_p_value_log10(mismatches_p),
        'mismatches_rbc': mismatches_rbc,

        'mapq_s': mapq_s,
        'mapq_p': handle_p_value_log10(mapq_p),
        'mapq_rbc': mapq_rbc,

        'mapq_mean': mapq_mean,
        'ref_mapq_mean':ref_mapq_mean,
        'alt_mapq_mean':alt_mapq_mean,

        'alt_UMI_avg_consistence_remove_single_read': alt_UMI_avg_consistence_remove_single_read,
        'alt_consistent_UMI_prop_relaxed_remove_single_read': alt_consistent_UMI_prop_relaxed_remove_single_read,
        'alt_consistent_UMI_prop_strict_remove_single_read': alt_consistent_UMI_prop_strict_remove_single_read,

        'read_number_s': read_number_s,
        'read_number_p': handle_p_value_log10(read_number_p),
        'read_number_rbc': read_number_rbc,

        'strand_bias_odds': strand_bias_odds,
        'strand_bias_p': handle_p_value_log10(strand_bias_p),

        'multi_mapper_odds': multi_mapper_odds,
        'multi_mapper_p': handle_p_value_log10(multi_mapper_p),
        'multi_map_prop':multi_map_prop,

        'softclip_length_s': softclip_length_s,
        'softclip_length_p': handle_p_value_log10(softclip_length_p),
        'softclip_length_rbc': softclip_length_rbc,

        'ref_allele_count': ref_allele_count,
        'consensus_ref_allele_count': consensus_ref_allele_count,
        'alt_allele_count': alt_allele_count,
        'consensus_alt_allele_count': consensus_alt_allele_count,
        'alt2_allele_count': alt2_allele_count,
        'consensus_alt2_allele_count': consensus_alt2_allele_count,

        'alt_vs_total_dp_paired_wilcoxon_rbc': alt_vs_total_dp_paired_wilcoxon_rbc,

        'UMI_end_remove_clip_s':UMI_end_remove_clip_s,
        'UMI_end_remove_clip_p':handle_p_value_log10(UMI_end_remove_clip_p),
        'UMI_end_remove_clip_rbc':UMI_end_remove_clip_rbc, 

        'reads_with_indel_odds': reads_with_indel_odds,
        'reads_with_indel_p': handle_p_value_log10(reads_with_indel_p),

        'ref_softclip_prop':ref_softclip_prop,
        'alt_softclip_prop':alt_softclip_prop,
        'softclip_prop':softclip_prop,

        'softclip_prop_odds':softclip_prop_odds,
        'softclip_prop_p':handle_p_value_log10(softclip_prop_p),

        'vaf_spot_mean':vaf_spot_mean,
        'vaf_spot_median':vaf_spot_median,
        'only_alt_mutant_spot_prop':only_alt_mutant_spot_prop, 

        'baseq_p_adj':handle_p_value_log10(baseq_p_adj),
        'ref_baseq1b_p_adj':handle_p_value_log10(ref_baseq1b_p_adj),
        'alt_baseq1b_p_adj':handle_p_value_log10(alt_baseq1b_p_adj),
        'querypos_p_adj':handle_p_value_log10(querypos_p_adj),
        'leftpos_p_adj':handle_p_value_log10(leftpos_p_adj),
        'seqpos_p_adj':handle_p_value_log10(seqpos_p_adj),
        'distance_to_end_p_adj':handle_p_value_log10(distance_to_end_p_adj),
        'UMI_end_p_adj':handle_p_value_log10(UMI_end_p_adj),
        'UMI_end_remove_clip_p_adj':handle_p_value_log10(UMI_end_remove_clip_p_adj),
        'mismatches_p_adj':handle_p_value_log10(mismatches_p_adj),
        'mapq_p_adj':handle_p_value_log10(mapq_p_adj),
        'read_number_p_adj':handle_p_value_log10(read_number_p_adj),
        'softclip_length_p_adj':handle_p_value_log10(softclip_length_p_adj)

    }
    
    # print(base_data['multi_mapper_odds'])
    if prior!="":
        fref="no"
        falt="no"
        prior_info=get_intergration_from_identifier_and_file(compare_pl_path,only_pos_identifier,prior,index_in_query="0 1 0 1")
        # print(prior_info)
        if prior_info!="":
            infos=prior_info.strip().split("\t")
            if "," not in ref:
                fref=str(infos[3+"ATCG".index(ref)])
            else:
                fref=infos[3+"ATCG".index(ref[0])] + "," + infos[3+"ATCG".index(ref[-1])]
            falt=infos[3+"ATCG".index(alt)]

        base_data['fref']=fref
        base_data['falt']=falt
    # features = features.reset_index()
    features_dict = features

    merged_data = {**{'identifier': identifier},**features_dict, **base_data}
    # print(merged_data['multi_mapper_odds'])
    return result_dict, merged_data



def main():
    # identifier_list=
    _global_features = None  # 占位符
    features=pd.read_csv(args.features,sep="\t")
    features.index=features['identifier']
    features_dict = features.to_dict(orient='index')  # 返回结构：{identifier: {col1: val1, ...}}
    def init_worker(features):
        global _global_features
        _global_features = features_dict

    # np.seterr(divide='raise')
    identifier_list=features.index.tolist()
    # out_file=os.path.join(outdir ,outname)
    import os
    # 假设 outdir 和 outname 已经定义
    outname = args.outname
    outdir = args.outdir
    if os.path.isabs(outname):
        out_file = os.path.join(outdir, os.path.basename(outname))
    else:
        out_file = os.path.join(outdir, outname)
    partial_func=partial(handel_identifier,args.bam,args.run_type,args.readLen,args.outdir,args.bins,args.prior)

    COLUMNS = []
    all_features = []

    with multiprocessing.Pool(processes=args.thread, initializer=init_worker, initargs=(features,)) as pool, open(out_file, "w") as f:
        for result in pool.imap(partial_func, identifier_list, chunksize=10):
            if result: 
                _, mutation_dict=result
                # print(mutation_dict)
                # all_features.append(result_dict)

                if mutation_dict is not None and mutation_dict!={}:
                    if COLUMNS!=[] or mutation_dict=={}:
                        pass
                    else:
                        COLUMNS=list(mutation_dict.keys())
                        pd.DataFrame(columns=COLUMNS).to_csv(f, header=True, index=False,sep='\t')

                    chunk_df = pd.DataFrame([[i for i in mutation_dict.values()]], columns=COLUMNS)
                    chunk_df.fillna("no").to_csv(f, header=False, mode="a", index=False,sep='\t')
            del result
            gc.collect()

    # with open(os.path.join(args.outdir, args.outname + ".pkl"), "wb") as pf:
    #     pickle.dump(all_features, pf)



# args.features="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/pre-classifier/scRNA/traning_data/features_manual/P4_normal.manual_addLabeling_addSampleid.identifier.feature.txt"
# thread=4
# outdir="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/pre-classifier/scRNA/traning_data/features_patched"
# outname="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/pre-classifier/scRNA/traning_data/features_random1w/P4_normal.random1w_patched.feature.txt"
# readLen=98
# bam="/storage/douyanmeiLab/yangzhirui/01.Data_download/06.skin/04.Analysis/01.cellranger_hg38/P4_normal/outs/possorted_genome_bam.bam"




## parameters
parser = argparse.ArgumentParser()
# parser.add_argument("--mutations","-m", required=False,default="",help="mutation identifier list")
parser.add_argument("--bins", required=False,default=100,type=int, help="only work for stereo-seq, if you want to combine the UMI in a bin level")
parser.add_argument("--features", required=True,help="feature file")
parser.add_argument("--prior", required=False,default="",help="prior file")
parser.add_argument("--thread", required=False,default=4,type=int, help="thread")
parser.add_argument("--outname", required=True,type=str, help="out file name")
parser.add_argument("--outdir", required=False,default="./",type=str, help="out dir name")
parser.add_argument("--type", dest='run_type',default="visium",choices=["visium","stereo","ST"],type=str, required=False, help="Your input sequence type")
parser.add_argument("--readLen", required=False,default=120,type=int, help="read length")
parser.add_argument("--bam",required=False,default="", help="bam_file")
args = parser.parse_args()
    
if __name__ == '__main__':
    main()  
