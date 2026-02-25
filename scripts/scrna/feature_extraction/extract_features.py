import shutil
import argparse
from collections import Counter, defaultdict
from datetime import datetime
from functools import partial
import multiprocessing
import os
import re
from statistics import mean
import statistics
import numpy as np
import pandas as pd
import subprocess
import uuid
import sys
from pybedtools import BedTool
import pysam
import scipy
import scipy.stats
from module.UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele, handle_cigar, handle_pos, handle_quality_matrix, handle_seq
from module.read_file import read_h5ad_file_2
from module.spatial_features import spatial_moran
from utils import calc_SegBias_for_one_sample, calc_vdb, combine_info_from_cigar, do_wilicox_sum_test, get_chr_size, get_indel_info, handle_p_value_log10, handle_posname,check_dir, judge_pos_in_indel, str_to_dict
from pyfaidx import Fasta


##input file: 
'''
file: phasing results
chr1    9271168 A       G       NA      chr1_9271218_G_T        H6PD    383     215     3       haplo>3 artifacts       A,G:125;A,T:87;G,G:1;G,T:2
chr1    9271195 A       G       NA      chr1_9271218_G_T        H6PD    372     269     4       haplo>3 artifacts       A,G:158;A,T:107;G,G:3;G,T:1
chr1    16618273        G       T       NA      chr1_16618314_G_C       CROCCP2 169     154     3       haplo>3 artifacts       G,G:82;G,C:69;T,G:2;T,C:1
chr1    16618274        C       A       NA      chr1_16618314_G_C       CROCCP2 170     155     4       haplo>3 artifacts       C,G:80;C,C:71;A,G:3;A,C:1

file: spatial feature.txt
#chr	pos	ref	alt	test_sig	early_mutation	late_mutation	verylate_mutation	ks_pvalue	
# moranI_pvalue	mutant_rate	mutant_rate_prob	mutant_rate_likelihood mutant_rate_vaf	mean_vaf	
# max_vaf	r_square	wilcoxon_pvalue	outlier_clusters	outlier_vaf	outlier_moranI_pvalue

'''


def refine_mean(in_list):
    try:
        return mean(in_list)
    except:
        print("wrong in refine mean for", in_list)
        return "NA"
    

def refine_median(in_list):
    try:
        return statistics.median(in_list)
    except:
        print("wrong in refine median for", in_list)
        return "NA"


def refine_diff(a,b):
    try:
        return a - b
    except:
        print("wrong in diff for: a:", a,"b:",b)
        return "NA"


class Features:
    def __init__(self,identifier):
        mut_chrom,mut_pos,mut_ref,mut_alt=handle_posname(identifier)
        #basic
        self.chrom = mut_chrom
        self.pos = mut_pos
        self.ref = mut_ref
        self.alt = mut_alt
        self.mut_origin = "NA"
        self.identifier=identifier
        # self.readLen=120

        #ind genotype
        self.spotNum="no"
        self.consensus_read_count="no"
        # self.refhom_likelihood="no"
        # self.althom_likelihood="no"
        # self.het_likelihood="no"
        self.mosaic_likelihood="no"
        self.AFind="no"

        #cluster info
        self.vaf_cluster_mean="no"
        self.vaf_cluster_std_dev="no"

        #spatial_feature
        self.r2="no"
        self.wilcoxon_p="no"
        self.wilcoxon_s="no"

        self.test_sig="no"
        self.KS_p="no"
        self.KS_s="no"
        self.MI_p="no"
        self.MI_s="no"
        self.mutant_rate="no"
        self.mutant_rate_prob="no"
        self.mutant_rate_likelihood="no"
        self.mutant_rate_vaf="no"
        self.mean_AFspot="no"
        self.max_AFspot="no"

        #spatial_feature: outlier cluter related
        self.outlier_moranI_pvalue="no"
        self.outlier_moranI_stat="no"

        self.outlier_clusters="no"
        self.outlier_vaf=0
        self.num_outlier_cluster=0

        #phase
        self.phase_condition="no"
        # self.phase_site=[]
        # self.haplotype=[]
        # self.phase_discordant=[]
        # self.phase_discordant_prop=[]
        # self.distancePhase=[]
        self.combine_nearest_phase_haplotype="no"
        self.combine_nearest_info_mutant_prop="no"
        self.combine_nearest_discordant_prop="no"
        self.combine_nearest_mut_origin="no"
        self.combine_nearest_phase_distance="no"

        self.no_combine_most_phase_haplotype="no"
        self.no_combine_most_info_mutant_prop="no"
        self.no_combine_most_discordant_prop="no"
        self.no_combine_most_mut_origin="no"
        self.no_combine_most_phase_distance="no"

        self.no_combine_nearest_phase_haplotype="no"
        self.no_combine_nearest_info_mutant_prop="no"
        self.no_combine_nearest_discordant_prop="no"
        self.no_combine_nearest_mut_origin="no"
        self.no_combine_nearest_phase_distance="no"

        self.combine_most_phase_haplotype="no"
        self.combine_most_info_mutant_prop="no"
        self.combine_most_discordant_prop="no"
        self.combine_most_mut_origin="no"
        self.combine_most_phase_distance="no"

        #gene_info
        self.gene_id = "no"
        self.gene_name = []
        self.gene_type = []
        self.transcript_name=[]
        self.exon_number=[]
        self.strand=[]
        self.distanceExon=[]
        # self.distanceExonIgnoreStrand=[]
        self.exonDesOrder=[]
        
        # read position
        self.ref_querypos_list="no"
        self.ref_querypos_num="no"
        self.alt_querypos_list="no"
        self.alt_querypos_num="no"

        self.left_ref_querypos_list_remove_clip="no"
        self.left_ref_querypos_num_remove_clip="no"
        self.left_alt_querypos_list_remove_clip="no"
        self.left_alt_querypos_num_remove_clip="no"

        self.right_ref_querypos_list_remove_clip="no"
        self.right_ref_querypos_num_remove_clip="no"
        self.right_alt_querypos_list_remove_clip="no"
        self.right_alt_querypos_num_remove_clip="no"
        
        #gene_count and exon_length
        self.exon_effect_len="no"
        self.total_exon_numbers="no"
        self.gene_count="no"
        self.fpkm="no"

        #expression variance
        self.MI_p_expression=[]
        self.MI_stat_expression=[]

        #mappability score
        self.mappability_score="no"

        #annotation from annovar
        self.anno="no"
        self.anno_gene="no"

        #fasta context and GCcontent
        # self.context="no"
        self.GCcontent="no"
        self.DNAMutationType="no"
        self.RNAMutationType="no"
        self.equal_to_previous_bases="no"
        self.cause_ploy_alt="no"

        #read-level info
        self.dp="no"
        self.alt_SpotNum="no"
        self.alt_dp="no"
        self.vaf="no"
        self.ref_hardclip_prop="no"
        self.alt_hardclip_prop="no"
        self.ref_softclip_prop="no"
        self.alt_softclip_prop="no"
        self.hardclip_length_s="no"
        self.hardclip_length_p="no"
        self.hardclip_prop_odds="no"
        self.hardclip_prop_p="no"
        self.softclip_length_s="no"
        self.softclip_length_p="no"
        self.softclip_prop_odds="no"
        self.softclip_prop_p="no"
        self.alt2_proportion="no"
        self.indel_proportion_SNPonly="no"
        self.mapq_difference="no"
        self.mapq_p="no"
        self.mapq_s="no"
        self.sb_p="no"
        self.sb_odds="no"
        self.mapper_odds="no"
        self.mapper_p="no"
        self.ref_multi_map_prop="no"
        self.alt_multi_map_prop="no"
        self.mismatches_p="no"
        self.mismatches_s="no"
        self.ref_mismatches_mean="no"
        self.alt_mismatches_mean="no"
        self.querypos_p="no"
        self.querypos_s="no"
        self.leftpos_p="no"
        self.leftpos_s="no"
        self.seqpos_p="no"
        self.seqpos_s="no"
        self.baseq_p="no"
        self.baseq_s="no"
        self.mean_distance_to_end="no"
        self.median_distance_to_end="no"
        self.mean_distance_to_end_remove_clip="no"
        self.median_distance_to_end_remove_clip="no"
        self.distance_to_end_p="no"
        self.distance_to_end_s="no"
        self.distance_to_end_remove_clip_p="no"
        self.distance_to_end_remove_clip_s="no"
        self.ref_baseq1b_p="no"
        self.ref_baseq1b_s="no"
        self.alt_baseq1b_p="no"
        self.alt_baseq1b_s="no"
        self.VDB="no"
        self.RPBZ="no"
        self.SGB="no"

        self.UMI_number_per_spot_s="no"
        self.UMI_number_per_spot_p="no"
        self.read_number_per_spot_s="no"
        self.read_number_per_spot_p="no"

        self.alt2_proportion_consensus="no"
        self.dp_consensus="no"
        self.alt_dp_consensus="no"
        self.vaf_consensus="no"

        self.reverse_dp="no"
        self.forward_dp="no"
        self.alt_reverse_dp="no"
        self.alt_forward_dp="no"
        self.major_read_strand="no"

        self.ref_ins_prop="no"
        self.ref_del_prop="no"
        self.ref_ins_major_prop="no"
        self.ref_del_major_prop="no"
        self.alt_ins_prop="no"
        self.alt_del_prop="no"
        self.alt_ins_major_prop="no"
        self.alt_del_major_prop="no"

        self.alt_UMI_consistence_prop="no"
        self.alt_consistence_hard_prop="no"
        self.alt_consistence_soft_prop="no"

        self.read_number_p="no"
        self.read_number_s="no"
        self.ref_read_number_perUMI_median="no"
        self.alt_read_number_perUMI_median="no"
        self.ref_read_number_perUMI_max="no"
        self.alt_read_number_perUMI_max="no"

        self.fref="no"
        self.falt="no"

        self.muts_in_cluster_p="no"
        self.muts_in_cluster_s="no"

    def add_info_from_ind_genotype(self, ind_geno_line):
        """
        #chrom	site	ID	germline	mutant	cluster	spot_number	consensus_read_count	genotype	p_mosaic	Gi	vaf

        """
        ind_geno_line=ind_geno_line.strip().split("\t")
        self.spotNum=ind_geno_line[6]
        self.consensus_read_count=ind_geno_line[7]
        if self.alt==ind_geno_line[4]:
            # self.refhom_likelihood=ind_geno_line[8]
            # self.althom_likelihood=ind_geno_line[9]
            # self.het_likelihood=ind_geno_line[10]
            self.mosaic_likelihood=ind_geno_line[9]
            self.AFind=ind_geno_line[11]

    def add_info_from_cluster_vaf(self,vaf_cluster_info):
        vaf_info=[]
        dps=[]
        mut_dps=[]
        for line in vaf_cluster_info.split("\n"):
            if line!="":
                sline=line.split("\t")
                vaf=sline[8]
                try:
                    vaf=float(vaf)
                    if vaf!=0.0:
                        vaf_info.append(vaf)
                except:
                    print(vaf)
                
                consensus_read_counts=sline[7].strip().split(",")
                dp=sum([int(i) for i in consensus_read_counts])
                alt_dp=consensus_read_counts["ATCG".index(self.alt)]

                dps.append(dp)
                mut_dps.append(alt_dp)

        mean = np.mean(vaf_info)
        std_dev = np.std(vaf_info, ddof=1)
        self.vaf_cluster_mean=mean
        self.vaf_cluster_std_dev=std_dev

        self.muts_in_cluster_s, self.muts_in_cluster_p=do_wilicox_sum_test(dps,mut_dps,type="list")


    def add_info_from_sf(self,sf_line):
        """
        get sf and AF info from spatial feature file:

        chr\tpos\tref\talt\ttest_sig\tearly_mutation\tlate_mutation\tverylate_mutation\t [0:7]
        ks_stat\tks_pvalue\t [8:9]
        moranI_stat\tmoranI_pvalue\t [10:11]
        mutant_rate\tmutant_rate_prob\tmutant_rate_likelihood\tmutant_rate_vaf\t [12:15]
        mean_vaf\tmax_vaf\t [16:17]
        r_squared\twilcoxon_stat\twilcoxon_pvalue\t [18:20]
        outlier_clusters\toutlier_vaf\toutlier_moranI_stat\toutlier_moranI_pvalue\n [21:24]

        """
        self.test_sig=sf_line[4]
        self.KS_s=sf_line[8]
        self.KS_p=sf_line[9]
        self.MI_s=sf_line[10]
        self.MI_p=sf_line[11]
        self.mutant_rate=sf_line[12]
        self.mutant_rate_prob=sf_line[13]
        self.mutant_rate_likelihood=sf_line[14]
        self.mutant_rate_vaf=sf_line[15]
        self.mean_AFspot=sf_line[16]
        self.max_AFspot=sf_line[17]

        self.r2=sf_line[18]
        self.wilcoxon_s=sf_line[19]
        self.wilcoxon_p=sf_line[20]
    
        self.outlier_clusters=sf_line[21]
        self.outlier_vaf=sf_line[22]
        self.outlier_moranI_stat=sf_line[23]
        self.outlier_moranI_pvalue=sf_line[24]
        self.num_outlier_cluster=len(self.outlier_clusters.strip().split(","))

        # self.outlier_vaf=sf_line[19]
        # self.nooutlier_vaf_ratio=sf_line[20]
        # self.nooutlier_prop_ratio=sf_line[21]
        # self.outlier_moranI_prob_pvalue=sf_line[22]
        # self.outlier_moranI_vaf_pvalue=sf_line[23]
        # self.num_outlier_cluster=len(self.outlier_clusters.strip().split(","))

    ####### this version the distance exon is list
    # def add_info_from_gff(self,sline):
    #     """
    #     gene_info, collected from gff file, download from gencode
    #     """
    #     try:
    #         exon_boundary_1=int(sline[6])-int(sline[1]); exon_boundary_2=int(sline[7])-int(sline[1])
    #         strand=sline[9]; Infos=sline[11]
    #         info_dict=handle_gff_info(Infos)
    #         exon_distance=exon_boundary_1 if abs(exon_boundary_1)<abs(exon_boundary_1) else exon_boundary_2
    #         if self.gene_id=="no":
    #             self.gene_id=[]

    #         if "transcript_id" in info_dict.keys():
    #             t_id=info_dict["transcript_id"].split(".")[0]
    #             print(self.transcript_name)
    #             print("t_id",t_id)
    #             if t_id in self.transcript_name:
    #                 self.gene_type.append(info_dict["gene_type"])
    #                 self.exon_number.append(info_dict["exon_number"])
    #                 self.strand.append(strand)
    #                 self.distanceExon.append(exon_distance)
    #     except:
    #         print(sline)

    def add_info_from_gff(self,sline):
        """
        gene_info, collected from gff file, download from gencode
        """

        try:
            exon_boundary_1=int(sline[6])-int(sline[1]); exon_boundary_2=int(sline[7])-int(sline[1])
            strand=sline[9]; Infos=sline[11]
            info_dict=handle_gff_info(Infos)
            exon_distance=abs(exon_boundary_1) if abs(exon_boundary_1)<abs(exon_boundary_1) else abs(exon_boundary_2)
            # print(exon_distance)
            # require the gene information from bam file
            if self.gene_id=="no":
                self.gene_id=[]

            if "transcript_id" in info_dict.keys() and strand==self.major_read_strand:
                # t_id=info_dict["transcript_id"].split(".")[0]
                self.gene_type.append(info_dict["gene_type"])
                # self.strand.append(strand)
                self.distanceExon.append(exon_distance)
                
            # self.distanceExonIgnoreStrand.append(exon_distance)
            
                # print(info_dict["transcript_id"],"====",self.transcript_name)
        except:
            print(sline)
            pass
    
    
    def add_info_from_wgEncodeGencodeExon(self,sline):
        """
        only distance exon
        #right:
        chr1    3073252 3073252 chr1    3073252 3074322 ENSMUST00000193812.1    AK016604.1      EMBL    ENSMUSE00001343744.1
        #no info
        chr22   22895504        22895504        .       -1      -1      .       .       .       .


        """
        try:
            if sline[3]!=".":
                exon_boundary_1=int(sline[4])-int(sline[1]); exon_boundary_2=int(sline[5])-int(sline[1])
                exon_distance=abs(exon_boundary_1) if abs(exon_boundary_1)<abs(exon_boundary_1) else abs(exon_boundary_2)
                self.distanceExon.append(exon_distance)
        except:
            print(sline)
            pass
        

    # def add_info_from_phase(self,phase_result,phase_type="combine"):
    #     """
    #     phase_info, from phase file
    #     chr1    9271143 A       G       NA      chr1_9271218_G_T        H6PD    73      52      6       haplo>3 artifacts       A,G:29;A,T:17;G,G:2;G,T:4
    #     chr1    9271168 A       G       NA      chr1_9271218_G_T        H6PD    81      64      3       haplo>3 artifacts       A,G:34;A,T:27;G,G:1;G,T:2
    #     chr1    9271195 A       G       NA      chr1_9271218_G_T        H6PD    88      78      4       haplo>3 artifacts       A,G:42;A,T:32;G,G:3;G,T:1
    #     """
    #     s_phase_result=phase_result.strip().split("\t")
    #     # print(s_phase_result)
    #     self.haplotype=s_phase_result[10]
    #     if self.haplotype=="haplo=3":
    #         details=s_phase_result[12]
    #         values = [int(item.split(":")[1]) for item in details.split(";")]
    #         pick_hSNP_value=values[0::2] if values[-2] < values[-1] else values[1::2] 
    #         # discordant=min(values[-2:]);discordant_prop=int(pick_hSNP_value[-1])/sum(pick_hSNP_value[0:-1])
    #         discordant=min(values[-2:]);discordant_prop=int((pick_hSNP_value[-1]))/int(sum(pick_hSNP_value[0:-1]))
    #         self.phase_discordant.append(discordant)
    #         self.phase_discordant_prop.append(discordant_prop)
    #     else:
    #         self.phase_discordant.append("no")
    #         self.phase_discordant_prop.append("no")
        
    #     phase_site=s_phase_result[5]
    #     _,phase_pos,_,_=handle_posname(phase_site)
    #     distancePhase=abs(int(phase_pos)-int(self.pos))
    #     self.distancePhase.append(distancePhase)

    #     if self.phase_condition=="no":
    #         self.phase_condition=[]
        
    #     self.phase_site=[]
    #     self.phase_condition.append(s_phase_result[10])
    #     self.phase_site.append(s_phase_result[5])
    #     self.mut_origin=s_phase_result[4]
    
    def add_info_from_phase(self,phase_result,phase_type="combine"):
        """
        phase_info, from phase file
        chr1    9271143 A       G       NA      chr1_9271218_G_T        H6PD    73      52      6       haplo>3 artifacts       A,G:29;A,T:17;G,G:2;G,T:4
        chr1    9271168 A       G       NA      chr1_9271218_G_T        H6PD    81      64      3       haplo>3 artifacts       A,G:34;A,T:27;G,G:1;G,T:2
        chr1    9271195 A       G       NA      chr1_9271218_G_T        H6PD    88      78      4       haplo>3 artifacts       A,G:42;A,T:32;G,G:3;G,T:1
        """
        s_phase_result=phase_result.strip().split("\n")
        s_phase_result=[line for line in s_phase_result if line !=""]
        if s_phase_result==[]:
            return
        distances=[int(line.strip().split("\t")[5].split("_")[1])-int(line.strip().split("\t")[1]) for line in s_phase_result]
        mutant_counts=[int(line.strip().split("\t")[9]) for line in s_phase_result]
        
        nearest_index=distances.index(min(distances))
        try:      
            most_mutant_index=mutant_counts.index(max(mutant_counts))
        except:
            print(s_phase_result)
        try:
            self.support_reads_prop_across_hSNPs=max(mutant_counts)/sum(mutant_counts)
        except:
            print(s_phase_result)
        for index,label in zip([nearest_index,most_mutant_index],["nearest","most"]):
            line=s_phase_result[index].strip().split("\t")
            haplotype=line[10]
            details=line[12]
            info_mutant_prop=int(line[9])/int(line[7]) # The mutant&infoSNP reads/allele&infSNP reads
            values = [int(item.split(":")[1]) for item in details.split(";")]
            pick_hSNP_value=values[0::2] if values[-2] < values[-1] else values[1::2] 
            try:
                discordant=min(values[-2:]);discordant_prop=int((pick_hSNP_value[-1]))/int(sum(pick_hSNP_value))
            except:
                discordant_prop=0;print(values,pick_hSNP_value)
            mut_origin=line[4]

            if label=="nearest" and phase_type=="combine":
                self.combine_nearest_phase_haplotype=haplotype
                self.combine_nearest_info_mutant_prop=info_mutant_prop
                self.combine_nearest_discordant_prop=discordant_prop
                self.combine_nearest_mut_origin=mut_origin
                self.combine_nearest_phase_distance=min(distances)
            elif label=="most" and phase_type=="no_combine":
                self.no_combine_most_phase_haplotype=haplotype
                self.no_combine_most_info_mutant_prop=info_mutant_prop
                self.no_combine_most_discordant_prop=discordant_prop
                self.no_combine_most_mut_origin=mut_origin
                self.no_combine_most_phase_distance=int(line[5].split("_")[1])-int(line[1])
            elif label=="nearest" and phase_type=="no_combine":
                self.no_combine_nearest_phase_haplotype=haplotype
                self.no_combine_nearest_info_mutant_prop=info_mutant_prop
                self.no_combine_nearest_discordant_prop=discordant_prop
                self.no_combine_nearest_mut_origin=mut_origin
                self.no_combine_nearest_phase_distance=min(distances)
            elif label=="most" and phase_type=="combine":
                self.combine_most_phase_haplotype=haplotype
                self.combine_most_info_mutant_prop=info_mutant_prop
                self.combine_most_discordant_prop=discordant_prop
                self.combine_most_mut_origin=mut_origin
                self.combine_most_phase_distance=int(line[5].split("_")[1])-int(line[1])

        if self.phase_condition=="no":
            self.phase_condition="yes" 

    def add_info_from_knownGene(self,knownGene_info):
        """
        exon_length,total exon numbers
        ENST00000456328.2       chr1    +       11868   14409   11868   11868   3       11868,12612,13220,      12227,12721,14409,              uc286dmu.1
        """
        s_transcript_info=knownGene_info.strip().split("\t")
        #print(s_transcript_info)
        total_exon_numbers=int(s_transcript_info[7])
        exon_start_info=[x.strip() for x in s_transcript_info[8].split(",") if x.strip()!='' ]
        exon_end_info=[x.strip() for x in s_transcript_info[9].split(",") if x.strip()!='' ]
        if len(exon_start_info) != len(exon_start_info):
            print(f"something wrong, the exon_strat_location is not equal to exon_end_location!\nStart:{exon_start_info}\nEnd:{exon_end_info}")
        
        exon_effect_len=0
        for exon_start, exon_end in zip(exon_start_info, exon_end_info):
            one_exon_length=int(exon_end)-int(exon_start)
            exon_effect_len+=one_exon_length

        if self.exon_effect_len=="no":
            self.exon_effect_len=[]
        if self.total_exon_numbers=="no":
            self.total_exon_numbers=[]

        self.exon_effect_len.append(exon_effect_len)
        self.total_exon_numbers.append(total_exon_numbers)

    def add_gene_count_info(self,gene_count_info):
        """
        ENSG00000000003.16      974
        """
        count=int(gene_count_info.strip().split("\t")[1])
        if self.gene_count=="no":
            self.gene_count=[]
        self.gene_count.append(count)


    def add_expression_info(self,total_gene_count):
        """
        cause the scRNA seq cannot extract the total length of mRNA, so the FPKM is not confident rather than gene count
        """
        if self.exon_effect_len!="no" and self.gene_count!="no":
            self.fpkm=[]
            for gene_count, exon_effect_len in zip(self.gene_count,self.exon_effect_len):
                fpkm = ((int(gene_count) / int(exon_effect_len)) / (total_gene_count / 1000))* 1000000
                # tpm = (gene_reads_count[i] / gene_length[i]) / sum([x / gene_length[j] for j, x in enumerate(gene_reads_count)]) * 1000000
                self.fpkm.append(fpkm)

    def add_mappability_info(self,mappability_info):
        """
        chr1	10157	10158	0.0
        """
        
        #if mappability_info!="" and len(mappability_info.strip().split("\n"))==1:
        if mappability_info!="":
            # print("1:",mappability_info)
            if type(self.mappability_score)==list:
                pass
            else:
                self.mappability_score=[]
            score=mappability_info.strip().split("\t")[3]
            self.mappability_score.append(score)
        # elif len(mappability_info.strip().split("\n"))>1:
        #     # print(">1:",mappability_info)
        #     # print("line:",mappability_info.strip())
        #     score=[]
        #     for line in mappability_info.strip().split("\n"):
        #         score.append(line.strip().split("\t")[3])
        #     self.mappability_score=",".join(score)

#     def add_annotation_from_annovar_table(self,annovar_info):
#         """
# chr1    2185380 2185380 T       A       exonic  FAAP20  .       nonsynonymous   SNV
# chr1    2310102 2310102 C       A       UTR3    SKI     NM_003036:c.*3337C>A    .       .
# chr1    3630128 3630128 G       A       downstream      TPRG1L;WRAP73   dist=642        .       .
# chr1    632663  632663  C       A       ncRNA_exonic    MIR12136        .       .       .
# chr1    633442  633442  C       A       upstream        MIR12136        dist=757        .       .
# chr1    633963  633963  T       G       intergenic      MIR12136;OR4F16 dist=1278;dist=51753    .       .
# chr1    19219002        19219002        G       A       ncRNA_intronic  EMC1-AS1        .       .       .
# chr1    153357876       153357876       G       T       UTR5    S100A9  NM_002965:c.-408G>T     .       .
# chr20   63735923        63735923        C       T       UTR5;UTR3       LIME1;ZGPAT     NM_001305654:c.-1627C>T;NM_001083113:c.*4C>T;NM_032527:c.*4C>T;NM_001195653:c.*4C>T;NM_001195654:c.*4C>T;NM_181485:c.*4C>T      .       .
#         """
#         # print(annovar_info)
#         s_annovar_info=annovar_info.strip().split("\t")
#         self.pos_anno=s_annovar_info[5]

#         if self.pos_anno=="exonic":
#             self.func_anno=s_annovar_info[8]


#     def add_annotation_from_annovar_all(self,annovar_info):
#         """
#         UTR5    RPS27(ENST00000651669.1:c.-35C>T)       chr1    153990762       153990762       C       T
#         UTR5    RPS27(ENST00000651669.1:c.-34C>T)       chr1    153990763       153990763       C       T
#         UTR3    RPS27(ENST00000651669.1:c.*54T>C)       chr1    153992147       153992147       T       C
#         UTR3    ARL8A(ENST00000272217.7:c.*1041T>A,ENST00000614750.1:c.*1090T>A),GPR37L1(ENST00000367282.6:c.*4870A>T)  chr1    202133426       202133426       A       T
#         intergenic      ENSG00000223649(dist=38643),SLC30A1(dist=40008) chr1    211531560       211531560       T       C
#         """
#         # print(annovar_info)
#         s_annovar_info=annovar_info.strip().split("\t")
#         self.pos_anno=s_annovar_info[0]


#     def add_annotation_from_annovar_exon(self,annovar_info):
#         """
#         line3   nonsynonymous SNV       RPS12:ENST00000230050.4:exon3:c.G73C:p.A25P,    chr6    132815030       132815030       G       C
#         """
#         s_annovar_info=annovar_info.strip().split("\t")
#         self.func_anno=s_annovar_info[1]

    def add_annotation_from_annovar(self,annovar_info):
        s_annovar_info=annovar_info.strip().split("\t")
        self.anno=s_annovar_info[1]
        self.anno_gene=s_annovar_info[2]



    def add_consensus_info(self,ind_count_info):
        """
        #chrom  pos     ID      ref     alt     cluster spot_num        umi_count       qA      qT      qC      qG
        GL000008.2      80303   .       C       G       bulk    2       0,0,1,1 NA      NA      48:1    43:1
        """
        s_ind_count_info=ind_count_info.strip().split("\t")
        # print(ind_count_info)
        qual_list=[]
        if "," in self.ref:
            for ref_allele in self.ref.split(","):
                ref_index="ATCG".index(ref_allele)
                qual_list.append(s_ind_count_info[8+ref_index])
            qual_str=",".join(qual_list)    
            qref=str_to_dict(qual_str)       
        else:
            ref_index="ATCG".index(self.ref)
            # print(ind_count_info,ref_index)
            qref=str_to_dict(s_ind_count_info[8+ref_index])
        alt_index="ATCG".index(self.alt)
        qalt=str_to_dict(s_ind_count_info[8+alt_index])
        self.baseq_p_consensus, self.baseq_s_consensus=do_wilicox_sum_test(qref,qalt,type="dict")
    
    def add_read_info(self,read_info_dict,dp):
        """
        read_info_dict[geno].keys:
        ["number_mismatch","is_reverse", "map_q", "is_indel", "baseq"
            "baseq1b", "rightpos_p", "leftpos_p","seqpos","querypos", "distance_to_end"
            "del_length", "del_distance", "del_num", "ins_distance", "ins_length","ins_num",
            "left_hardclip","right_hardclip","hardclip_length", "left_softclip", "right_softclip", "softclip_length"]    
        """
        def simple_to_get_list(read_info_dict,ref_allele,alt_allele,var):
            # print(ref_allele,alt_allele, ref_allele.split(","))
            return_ref_list=[int(k) for allele in ref_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
            return_alt_list=[int(k) for allele in alt_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
            return return_ref_list,return_alt_list
            
        get_list=partial(simple_to_get_list,read_info_dict,self.ref,self.alt)


        # gene id and trancript id
        GeneID_list=[]
        GeneName_list=[]
        TransID_list=[]
        for allele in "ATCG":
            GeneID_list += read_info_dict[allele]["GeneID_list"]
            GeneName_list += read_info_dict[allele]["GeneName_list"]
            TransID_list += read_info_dict[allele]["TransID_list"]

        def tidy_result(in_list):
            a_list=list(set(in_list))
            if "no" in a_list:
                a_list.remove("no")
            if a_list is None:
                a_list=[]
            b_list=[]    
            for m in a_list:
                if ";" in m or "," in m:
                    m_a=re.split(";|,",m)
                    b_list.append(m_a)
                else:
                    b_list.append(m)
            a_list=b_list
            return a_list
        self.gene_id=tidy_result(GeneID_list)
        self.gene_name=tidy_result(GeneName_list)
        self.transcript_name=tidy_result(TransID_list)

        ref_mismatches, alt_mismatches=get_list("number_mismatch")
        refine_alt_mismatches=[int(i)-1 if i!=0 else 0 for i in alt_mismatches]
        
        #number_multi_mapper
        ref_mappers, alt_mappers=get_list("number_mapper")
        refine_ref_mappers_uniq = len([i for i in ref_mappers if i == 1]); refine_ref_mappers_multi=len(ref_mappers)-refine_ref_mappers_uniq
        refine_alt_mappers_uniq = len([i for i in alt_mappers if i == 1]); refine_alt_mappers_multi=len(alt_mappers)-refine_alt_mappers_uniq
        self.mapper_odds, self.mapper_p = scipy.stats.fisher_exact([[refine_alt_mappers_multi, refine_alt_mappers_uniq ],[refine_ref_mappers_multi, refine_ref_mappers_uniq]])
        self.ref_multi_map_prop=refine_ref_mappers_multi/len(ref_mappers) if len(ref_mappers)!=0 else "NA"
        self.alt_multi_map_prop=refine_alt_mappers_multi/len(alt_mappers) if len(alt_mappers)!=0 else "NA"

        self.ref_mismatches_mean=refine_mean(ref_mismatches)
        self.alt_mismatches_mean=refine_mean(refine_alt_mismatches)
        self.mismatches_s,self.mismatches_p=do_wilicox_sum_test(ref_mismatches,refine_alt_mismatches,type="list")
        
        ref_is_reverse, alt_is_reverse=get_list("is_reverse")
        refine_ref_is_reverse = len([i for i in ref_is_reverse if i == 1]); refine_ref_is_forward = len(ref_is_reverse)-refine_ref_is_reverse
        refine_alt_is_reverse = len([i for i in alt_is_reverse if i == 1]); refine_alt_is_forward = len(alt_is_reverse)-refine_alt_is_reverse
        self.sb_odds, self.sb_p = scipy.stats.fisher_exact([[refine_alt_is_reverse, refine_alt_is_forward ],[refine_ref_is_reverse, refine_ref_is_forward]])

        ref_map_q, alt_map_q=get_list("map_q")
        self.mapq_s,self.mapq_p=do_wilicox_sum_test(ref_map_q,alt_map_q,method="greater",type="list")
        self.mapq_difference=refine_diff(refine_mean(ref_map_q),refine_mean(alt_map_q))

        self.dp=dp
        alt_dp= int(read_info_dict[self.alt]["dp"])
        self.alt_dp=alt_dp
        if dp!=0:
            self.vaf=int(alt_dp)/int(dp)
        else:
            self.vaf=0
        if dp!=0:
            self.indel_proportion_SNPonly=sum(read_info_dict["del"]["is_indel"])/dp
        else:
            self.indel_proportion_SNPonly=0

        self.reverse_dp=sum([int(read_info_dict[allele]["reverse_dp"]) for allele in "ATCG"])
        self.forward_dp=sum([int(read_info_dict[allele]["forward_dp"]) for allele in "ATCG"])
        self.alt_reverse_dp=int(read_info_dict[self.alt]["reverse_dp"])
        self.alt_forward_dp=int(read_info_dict[self.alt]["forward_dp"])

        if self.reverse_dp>=self.forward_dp:
            self.major_read_strand="-"
        elif self.reverse_dp<self.forward_dp:
            self.major_read_strand="+"

        ref_baseq, alt_baseq=get_list("baseq")
        self.baseq_s,self.baseq_p=do_wilicox_sum_test(ref_baseq,alt_baseq,method="greater",type="list")

        ref_baseq1b,alt_baseq1b=get_list("baseq1b")
        self.ref_baseq1b_s, self.ref_baseq1b_p=do_wilicox_sum_test(ref_baseq,ref_baseq1b,method="greater",type="list")
        self.alt_baseq1b_s, self.alt_baseq1b_p=do_wilicox_sum_test(alt_baseq,alt_baseq1b,method="greater",type="list")
        # print(ref_baseq,ref_baseq1b)
        if self.dp!=0:
            self.alt2_proportion=int(max([read_info_dict[allele]["dp"] for allele in "ATCG" if allele not in self.ref+self.alt]))/self.dp

        # ref_rightpos,alt_rightpos=get_list("rightpos_p")
        # _,self.rightpos_p=do_wilicox_sum_test(ref_rightpos,alt_rightpos,type="list")

        ref_leftpos,alt_leftpos=get_list("leftpos_p")
        self.leftpos_s,self.leftpos_p=do_wilicox_sum_test(ref_leftpos,alt_leftpos,type="list")

        ref_seqpos,alt_seqpos=get_list("seqpos")
        self.seqpos_s,self.seqpos_p=do_wilicox_sum_test(ref_seqpos,alt_seqpos,type="list")

        ref_querypos,alt_querypos=get_list("querypos")
        self.querypos_s,self.querypos_p=do_wilicox_sum_test(ref_querypos,alt_querypos,type="list")
        self.ref_querypos_list = ",".join([str(i) for i in ref_seqpos]); self.ref_querypos_num=len(set(ref_seqpos))
        self.alt_querypos_list = ",".join([str(i) for i in alt_seqpos]); self.alt_querypos_num=len(set(alt_seqpos))
     
        left_ref_edist,left_alt_edist=get_list("left_read_edist")
        # print("left_ref_edist\n",set(left_ref_edist))
        # print("left_alt_edist\n",set(left_alt_edist))
        self.left_ref_querypos_list_remove_clip = ",".join([str(i) for i in left_ref_edist])
        self.left_ref_querypos_num_remove_clip=len(set(left_ref_edist))
        self.left_alt_querypos_list_remove_clip = ",".join([str(i) for i in left_alt_edist])
        self.left_alt_querypos_num_remove_clip=len(set(left_alt_edist))

        right_ref_edist,right_alt_edist=get_list("right_read_edist")
        # print("right_ref_edist\n",set(right_ref_edist))
        # print("right_alt_edist\n",set(right_alt_edist))
        self.right_ref_querypos_list_remove_clip = ",".join([str(i) for i in right_ref_edist])
        self.right_ref_querypos_num_remove_clip=len(set(right_ref_edist))
        self.right_alt_querypos_list_remove_clip = ",".join([str(i) for i in right_alt_edist])
        self.right_alt_querypos_num_remove_clip=len(set(right_alt_edist))

        # self.queryposNum_odds, self.queryposNum_p = scipy.stats.fisher_exact([[refine_alt_mappers_multi, refine_alt_mappers_uniq ],[refine_ref_mappers_multi, refine_ref_mappers_uniq]])

        distance_list=[]
        distance_list_remove_clip=[]
        for allele in "ATCG":
            distance_list+=read_info_dict[allele]["distance_to_end"]
            distance_list_remove_clip+=read_info_dict[allele]["distance_to_end_remove_clip"]
        self.mean_distance_to_end=refine_mean(distance_list)
        self.median_distance_to_end=refine_median(distance_list)
        self.mean_distance_to_end_remove_clip=refine_mean(distance_list_remove_clip)
        self.median_distance_to_end_remove_clip=refine_median(distance_list_remove_clip)

        ref_distance_to_end,alt_distance_to_end=get_list("distance_to_end")
        self.distance_to_end_s,self.distance_to_end_p=do_wilicox_sum_test(ref_distance_to_end,alt_distance_to_end,type="list")

        ref_distance_to_end_remove_clip,alt_distance_to_end_remove_clip=get_list("distance_to_end_remove_clip")
        self.distance_to_end_remove_clip_s,self.distance_to_end_remove_clip_p=do_wilicox_sum_test(ref_distance_to_end_remove_clip,alt_distance_to_end_remove_clip,type="list")

        ## NOTE: save softclip_length and plot them
        ref_hardclip_length,alt_hardclip_length=get_list("hardclip_length")
        self.hardclip_length_s,self.hardclip_length_p=do_wilicox_sum_test(ref_hardclip_length,alt_hardclip_length,type="list")
        try:
            self.ref_hardclip_prop=len([x for x in ref_hardclip_length if x > 10])/len(ref_hardclip_length)
            self.alt_hardclip_prop=len([x for x in alt_hardclip_length if x > 10])/len(alt_hardclip_length)
            # self.ref_hardclip_length,self.alt_hardclip_length=ref_hardclip_length,alt_hardclip_length
            ref_hard_count=len([x for x in ref_hardclip_length if x > 10]); ref_no_hard_count=len(ref_hardclip_length)-ref_hard_count
            alt_hard_count=len([x for x in alt_hardclip_length if x > 10]); alt_no_hard_count=len(alt_hardclip_length)-alt_hard_count
            self.hardclip_prop_odds, self.hardclip_prop_p=scipy.stats.fisher_exact([[alt_hard_count, alt_no_hard_count ],[ref_hard_count, ref_no_hard_count]])
    
        except:
            print(self.identifier, "dose not have hard_clip_info")

        ref_softclip_length,alt_softclip_length=get_list("softclip_length")
        self.softclip_length_s,self.softclip_length_p=do_wilicox_sum_test(ref_softclip_length,alt_softclip_length,type="list")
        try:
            self.ref_softclip_prop=len([x for x in ref_softclip_length if x > 10])/len(ref_softclip_length)
            self.alt_softclip_prop=len([x for x in alt_softclip_length if x > 10])/len(alt_softclip_length)
            # self.ref_softclip_length,self.alt_softclip_length=ref_softclip_length,alt_softclip_length
            ref_soft_count=len([x for x in ref_softclip_length if x > 10]); ref_no_soft_count=len(ref_softclip_length)-ref_soft_count
            alt_soft_count=len([x for x in alt_softclip_length if x > 10]); alt_no_soft_count=len(alt_softclip_length)-alt_soft_count
            self.softclip_prop_odds, self.softclip_prop_p=scipy.stats.fisher_exact([[alt_soft_count, alt_no_soft_count ],[ref_soft_count, ref_no_soft_count]])
            
        except:
            print(self.identifier, "dose not have soft_clip_info")
        # def simple_to_get_list_float(read_info_dict,ref_allele,alt_allele,var):
        #     # print(ref_allele,alt_allele, ref_allele.split(","))
        #     return_ref_list=[k for allele in ref_allele.split(",") for k in read_info_dict[allele][var] if k !=0.0]
        #     return_alt_list=[k for allele in alt_allele.split(",") for k in read_info_dict[allele][var] if k !=0.0]
        #     return return_ref_list,return_alt_list

        # get_list_float=partial(simple_to_get_list_float,read_info_dict,self.ref,self.alt)
        #UMI_consistence_prop
        # _,alt_UMI_consistence=get_list_float("UMI_consistence_prop")
        alt_UMI_consistence=0; total_UMI_contain_alt=0; alt_consistence_hard=0; alt_consistence_soft=0
        for allele in self.alt.split(","):
            for k in read_info_dict[allele]["UMI_consistence_prop"]:
                if k !=0.0:
                    alt_UMI_consistence+=k
                    total_UMI_contain_alt+=1
                    if k==1.0:
                        alt_consistence_hard+=1
                        alt_consistence_soft+=1
                    elif k>=0.75:
                        alt_consistence_soft+=1
        try:
            alt_UMI_consistence_prop=alt_UMI_consistence/total_UMI_contain_alt
            alt_consistence_hard_prop=alt_consistence_hard/total_UMI_contain_alt
            alt_consistence_soft_prop=alt_consistence_soft/total_UMI_contain_alt
            self.alt_UMI_consistence_prop=alt_UMI_consistence_prop
            self.alt_consistence_hard_prop=alt_consistence_hard_prop
            self.alt_consistence_soft_prop=alt_consistence_soft_prop
        except:
            print(self.identifier,"is wrong in consistence.")
        try:
            ref_read_number_perUMI,alt_read_number_perUMI=get_list("read_number_per_UMI")
            self.read_number_s,self.read_number_p=do_wilicox_sum_test(ref_read_number_perUMI,alt_read_number_perUMI,type="list")
            self.ref_read_number_perUMI_median=refine_median(ref_read_number_perUMI)
            self.alt_read_number_perUMI_median=refine_median(alt_read_number_perUMI)
            self.ref_read_number_perUMI_max=max(ref_read_number_perUMI)
            self.alt_read_number_perUMI_max=max(alt_read_number_perUMI)
        except:
            print(self.identifier,"does not have perUMI info")
        
        try:
            ref_read_number_per_spot,alt_read_number_per_spot=get_list("read_number_per_spot")
            self.read_number_per_spot_s,self.read_number_per_spot_p = do_wilicox_sum_test(ref_read_number_per_spot, alt_read_number_per_spot,method="two-sided",type="list")
        except:
            print(self.identifier,"does not have UMI_number_per_spot info")
        
        # alt_consistence_prop=sum(alt_UMI_consistence)/sum([1 for i in read_info_dict[self.alt]["UMI_consistence_prop"] if i !=0.0])
        # total_UMI_contain_alt=[i for i in read_info_dict[self.alt]["UMI_consistence_prop"] if i !=0.0]
        # alt_consistence_hard=sum([1 for i in total_UMI_contain_alt if i ==1.0])/len(total_UMI_contain_alt)
        # alt_consistence_soft=sum([1 for i in total_UMI_contain_alt if i >=0.75])/len(total_UMI_contain_alt)       
        
        # "del_length", "del_distance", "del_num", "ins_distance", "ins_length","ins_num",

        def simple_to_get_list_facing_list(read_info_dict,ref_allele,alt_allele,var):
            # print(read_info_dict[ref_allele][var])
            return_ref_list=[i for allele in ref_allele.split(",") for k in read_info_dict[allele][var] for i in k]
            return_alt_list=[i for allele in alt_allele.split(",") for k in read_info_dict[allele][var] for i in k]
            # return_ref_list=[str(k) for allele in ref_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
            # return_alt_list=[str(k) for allele in alt_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
            return return_ref_list,return_alt_list
            
        get_list_face_list=partial(simple_to_get_list_facing_list,read_info_dict,self.ref,self.alt)
        #ref_ins_distance: [[],["no",1],[2]]
        ref_ins_num,alt_ins_num=get_list("ins_num")
        ref_ins_length,alt_ins_length=get_list_face_list("ins_length")
        ref_ins_distance,alt_ins_distance=get_list_face_list("ins_distance")
        self.ref_ins_num,self.alt_ins_num=ref_ins_num,alt_ins_num
        self.ref_ins_length,self.alt_ins_length=ref_ins_length,alt_ins_length
        self.ref_ins_distance,self.alt_ins_distance=ref_ins_distance,alt_ins_distance

        ref_del_num,alt_del_num=get_list("del_num")
        ref_del_length,alt_del_length=get_list_face_list("del_length")
        ref_del_distance,alt_del_distance=get_list_face_list("del_distance")
        self.ref_del_num,self.alt_del_num=ref_del_num,alt_del_num
        self.ref_del_length,self.alt_del_length=ref_del_length,alt_del_length
        self.ref_del_distance,self.alt_del_distance=ref_del_distance,alt_del_distance

        def get_sum_except_0(indel_number_list):
            count_dict=Counter(indel_number_list)
            sum_value=sum([int(count_dict[key]) for key in count_dict.keys() if key not in ["0",0]])
            total=sum(count_dict.values())
            try:
                return sum_value/total
            except:
                return 0
        
        def get_major_prop_except_no(ref_del_length):
            count_dict=Counter(ref_del_length);del count_dict["no"]
            if count_dict!={}:
                most_common_key = count_dict.most_common(1)[0][0]
                most_common_count=count_dict[most_common_key]
                indel_num_except_no=sum(count_dict.values())
                most_common_prop=int(most_common_count)/int(indel_num_except_no)
                return most_common_key,most_common_prop
            else:
                return "",0
        
        sum_ref_ins_num=get_sum_except_0(self.ref_ins_num)
        sum_ref_del_num=get_sum_except_0(self.ref_del_num)
        sum_alt_ins_num=get_sum_except_0(self.alt_ins_num)
        sum_alt_del_num=get_sum_except_0(self.alt_del_num)

        self.ref_ins_prop=sum_ref_ins_num
        self.ref_del_prop=sum_ref_del_num
        _,self.ref_ins_major_prop=get_major_prop_except_no(self.ref_ins_distance)
        _,self.ref_del_major_prop=get_major_prop_except_no(self.ref_del_distance)
        self.alt_ins_prop=sum_alt_ins_num
        self.alt_del_prop=sum_alt_del_num
        _,self.alt_ins_major_prop=get_major_prop_except_no(self.alt_ins_distance)
        _,self.alt_del_major_prop=get_major_prop_except_no(self.alt_del_distance)

        # features from bcftools: VDB   
        self.VDB=calc_vdb(read_info_dict[self.alt]["edist"],self.readLen,self.readLen)
        ref_pos,alt_pos=get_list("epos")
        self.RPBZ,_=do_wilicox_sum_test(ref_pos,alt_pos,method="two-sided",type="list")
        # self.MQBZ,_=calc_mwu_biasZ(bca->ref_mq,  bca->alt_mq, bca->nqual,1,1) #mapQ have been done before
        # BQBZ= calc_mwu_biasZ(bca->ref_bq,  bca->alt_bq, bca->nqual,0,1) #baseQ have been done before
        # MQSBZ= calc_mwu_biasZ(bca->fwd_mqs, bca->rev_mqs, bca->nqual,0,1) # this have no sense, cause all RNA-seq reads are same direction
        ref_dp_list=[read_info_dict[allele]["dp"] for allele in self.ref.split(",")]
        alt_dp_list=[read_info_dict[allele]["dp"] for allele in self.alt.split(",")]
        self.SGB=calc_SegBias_for_one_sample(sum(ref_dp_list),sum(alt_dp_list))

        self.alt_SpotNum=read_info_dict[self.alt]["GenoSpotNum"]


    def add_context_info(self,GCcontent,DNAMutationType,RNAMutationType, equal_to_previous_bases, cause_ploy_alt):
        # self.context=context
        self.GCcontent=GCcontent
        self.DNAMutationType=DNAMutationType
        self.RNAMutationType=RNAMutationType
        self.equal_to_previous_bases=equal_to_previous_bases
        self.cause_ploy_alt=cause_ploy_alt

    def add_info_from_barcode_dir(self, ref_spot_dp_list, alt_spot_dp_list):
        self.UMI_number_per_spot_s,self.UMI_number_per_spot_p = do_wilicox_sum_test(ref_spot_dp_list, alt_spot_dp_list,method="two-sided",type="list")

    def add_prior(self, prior_info):
        infos=prior_info.strip().split("\t")
        print(infos)
        # print("ATCG".index(self.ref))
        if "," not in self.ref:
            self.fref=str(infos[3+"ATCG".index(self.ref)])
        else:
            self.fref=infos[3+"ATCG".index(self.ref[0])] + "," + infos[3+"ATCG".index(self.ref[-1])]
        self.falt=infos[3+"ATCG".index(self.alt)]
        print(self.fref,self.falt)

    def expand_features(self):
        """
        after combine the features from different files, some features can be expanded
        such as: 
        exonDesOrder: sort descending of exon location
        
        """
        if self.total_exon_numbers!="no" and self.exon_number!="no":
            if self.exonDesOrder=="no":
                self.exonDesOrder=[]
            for total_exon_numbers,exon_number in zip(self.total_exon_numbers,self.exon_number):
                index_num=int(total_exon_numbers)-int(exon_number)+1
                self.exonDesOrder.append(index_num)

        if self.consensus_read_count!="no" and self.ref!="no" and self.alt!="no":
            other_count_dict={base:int(count) for base,count in zip("ATCG",self.consensus_read_count.split(",")) if base not in self.ref+self.alt}
            # print(self.consensus_read_count,"===", other_count_dict)
            alt2_count_consensus=max(list(other_count_dict.values()),default=0)
            consensus_dp=sum([int(i) for i in self.consensus_read_count.split(",")])
            self.alt2_proportion_consensus=alt2_count_consensus/consensus_dp
            self.dp_consensus=consensus_dp
            self.alt_dp_consensus=self.consensus_read_count.split(",")["ATCG".index(self.alt)]
            self.vaf_consensus=int(self.alt_dp_consensus)/self.dp_consensus

        if self.strand!=[]:
            self.strand=list(set(self.strand))

        if self.gene_type!=[]:
            self.gene_type=list(set(self.gene_type))
    
    def to_dict(self):
        # print(self.MI_mutant_p,type(self.MI_mutant_p),self.outlier_moranI_prob_pvalue,type(self.outlier_moranI_prob_pvalue),self.outlier_moranI_vaf_pvalue,type(self.outlier_moranI_vaf_pvalue))
        back_dict= {
            'identifier': self.identifier,
            # read level
            'querypos_p': handle_p_value_log10(self.querypos_p),
            'querypos_s': self.querypos_s,
            'leftpos_p': handle_p_value_log10(self.leftpos_p),
            'leftpos_s': self.leftpos_s,
            'seqpos_p': handle_p_value_log10(self.seqpos_p),
            'seqpos_s': self.seqpos_s,
            'mean_distance_to_end': self.mean_distance_to_end,
            'median_distance_to_end':self.median_distance_to_end,
            'mean_distance_to_end_remove_clip':self.mean_distance_to_end_remove_clip,
            'median_distance_to_end_remove_clip':self.median_distance_to_end_remove_clip,
            'distance_to_end_p':handle_p_value_log10(self.distance_to_end_p),
            'distance_to_end_s':self.distance_to_end_s,
            'distance_to_end_remove_clip_p':handle_p_value_log10(self.distance_to_end_remove_clip_p),
            'distance_to_end_remove_clip_s':self.distance_to_end_remove_clip_s,
            # 'querypos_p': (self.querypos_p),
            # 'leftpos_p': (self.leftpos_p),
            # 'seqpos_p': (self.seqpos_p),
            # 'distanceBoundary': self.distanceBoundary,

            'baseq_p': handle_p_value_log10(self.baseq_p),
            'baseq_s': self.baseq_s,
            'ref_baseq1b_p': handle_p_value_log10(self.ref_baseq1b_p),
            'ref_baseq1b_s': self.ref_baseq1b_s,
            'alt_baseq1b_p': handle_p_value_log10(self.alt_baseq1b_p),
            'alt_baseq1b_s': self.alt_baseq1b_s,
            # 'context': self.context,
            'DNAMutationType': self.DNAMutationType,
            'RNAMutationType': self.RNAMutationType,
            'equal_to_previous_bases': self.equal_to_previous_bases,
            'cause_ploy_alt': self.cause_ploy_alt,
            'mapper_odds': handle_p_value_log10(self.mapper_odds),
            'mapper_p': handle_p_value_log10(self.mapper_p),
            'ref_multi_map_prop': self.ref_multi_map_prop,
            'alt_multi_map_prop': self.alt_multi_map_prop,
            'ref_mismatches_mean': self.ref_mismatches_mean,
            'alt_mismatches_mean': self.alt_mismatches_mean,
            'mismatches_p': handle_p_value_log10(self.mismatches_p),
            'mismatches_s': self.mismatches_s,
            'sb_p': handle_p_value_log10(self.sb_p),
            'sb_odds': handle_p_value_log10(self.sb_odds),
            'mapq_p': handle_p_value_log10(self.mapq_p),
            'mapq_s': self.mapq_s,
            'mapq_difference': self.mapq_difference,
            'dp': self.dp,
            'alt_SpotNum': self.alt_SpotNum,
            'alt_dp': self.alt_dp,
            'vaf': self.vaf,
            'dp_consensus':self.dp_consensus,
            'alt_dp_consensus':self.alt_dp_consensus,
            'vaf_consensus': self.vaf_consensus,
            'VDB': self.VDB,
            'RPBZ':self.RPBZ,
            'SGB':self.SGB,
            # 'dpUMI_p': self.dpUMI_p,
            # 'dp_diff': self.feadp_diffture,
            "reverse_dp": self.reverse_dp,
            "forward_dp": self.forward_dp,
            "alt_reverse_dp": self.alt_reverse_dp,
            "alt_forward_dp": self.alt_forward_dp,
            "major_read_strand": self.major_read_strand,
            'ref_hardclip_prop': self.ref_hardclip_prop,
            'alt_hardclip_prop': self.alt_hardclip_prop,
            'ref_softclip_prop': self.ref_softclip_prop,
            'alt_softclip_prop': self.alt_softclip_prop,
            'hardclip_length_s': self.hardclip_length_s,
            'hardclip_length_p': handle_p_value_log10(self.hardclip_length_p),
            'hardclip_prop_odds': self.hardclip_prop_odds,
            'hardclip_prop_p': self.hardclip_prop_p,
            'softclip_length_s': self.softclip_length_s,
            'softclip_length_p': handle_p_value_log10(self.softclip_length_p),
            'softclip_prop_odds': self.softclip_prop_odds,
            'softclip_prop_p': handle_p_value_log10(self.softclip_prop_p),
            'indel_proportion_SNPonly': self.indel_proportion_SNPonly,
            'alt2_proportion': self.alt2_proportion,
            'alt2_proportion_consensus': self.alt2_proportion_consensus,
            'UMI_number_per_spot_s': self.UMI_number_per_spot_s,
            'UMI_number_per_spot_p': handle_p_value_log10(self.UMI_number_per_spot_p),
            'read_number_per_spot_s': self.read_number_per_spot_s,
            'read_number_per_spot_p': handle_p_value_log10(self.read_number_per_spot_p),
            # 'nearSNV': self.nearSNV,
            # 'nearIND_DEL': self.nearIND_DEL,
            # 'nearIND_DEL_length': self.nearIND_DEL_length,
            # 'nearIND_DEL_distance': self.nearIND_DEL_distance,
            # 'splicing': self.splicing,
            # UMI consistence:
            'alt_UMI_consistence_prop': self.alt_UMI_consistence_prop,
            'alt_consistence_hard_prop': self.alt_consistence_hard_prop,
            'alt_consistence_soft_prop': self.alt_consistence_soft_prop,
            'read_number_p': handle_p_value_log10(self.read_number_p),
            'read_number_s':self.read_number_s,
            'ref_read_number_perUMI_median':self.ref_read_number_perUMI_median,
            'alt_read_number_perUMI_median':self.alt_read_number_perUMI_median,
            'ref_read_number_perUMI_max':self.ref_read_number_perUMI_max,
            'alt_read_number_perUMI_max':self.alt_read_number_perUMI_max,

            # site info
            'spotNum': self.spotNum,
            'combine_nearest_phase_haplotype':self.combine_nearest_phase_haplotype,
            'combine_nearest_info_mutant_prop':self.combine_nearest_info_mutant_prop,
            'combine_nearest_discordant_prop':self.combine_nearest_discordant_prop,
            'combine_nearest_phase_distance':self.combine_nearest_phase_distance,

            'no_combine_most_phase_haplotype':self.no_combine_most_phase_haplotype,
            'no_combine_most_info_mutant_prop':self.no_combine_most_info_mutant_prop,
            'no_combine_most_discordant_prop':self.no_combine_most_discordant_prop,
            'no_combine_most_phase_distance':self.no_combine_most_phase_distance,

            'no_combine_nearest_phase_haplotype':self.no_combine_nearest_phase_haplotype,
            'no_combine_nearest_info_mutant_prop':self.no_combine_nearest_info_mutant_prop,
            'no_combine_nearest_discordant_prop':self.no_combine_nearest_discordant_prop,
            'no_combine_nearest_phase_distance':self.no_combine_nearest_phase_distance,

            'combine_most_phase_haplotype':self.combine_most_phase_haplotype,
            'combine_most_info_mutant_prop':self.combine_most_info_mutant_prop,
            'combine_most_discordant_prop':self.combine_most_discordant_prop,
            'combine_most_phase_distance':self.combine_most_phase_distance,

            'mappabilityScore': self.mappability_score,
            'GCcontent': self.GCcontent,
            'gene_type': self.gene_type,
            'anno_gene': self.anno_gene ,
            'anno': self.anno,
            # 'func_anno': self.func_anno,
            # 'RNAediting': self.RNAediting,
            # 'GeneNumber': self.GeneNumber,
            'geneStrand': self.major_read_strand,
            'distanceExon': self.distanceExon,
            # 'distanceExonIgnoreStrand':self.distanceExonIgnoreStrand,
            'AFind': self.AFind,
            'mean_AFspot': self.mean_AFspot,
            'max_AFspot': self.max_AFspot,

            'vaf_cluster_mean':self.vaf_cluster_mean,
            'vaf_cluster_std_dev':self.vaf_cluster_std_dev,

            # gene expression
            'MI_p_expression': handle_p_value_log10(self.MI_p_expression),
            'MI_s_expression': self.MI_stat_expression,

            # statistic info
            'mosaic_likelihood': self.mosaic_likelihood,
            # 'het_likelihood': self.het_likelihood,
            # 'althom_likelihood': self.althom_likelihood,
            # 'refhom_likelihood': self.refhom_likelihood,

            # spatial related
            'r2': self.r2,
            'wilcoxon_p':handle_p_value_log10(self.wilcoxon_p),
            'wilcoxon_s':self.wilcoxon_s,
            'KS_p': handle_p_value_log10(self.KS_p),
            'KS_s': self.KS_s,
            'MI_p': handle_p_value_log10(self.MI_p),
            'MI_s': self.MI_s,
            'sf_test_sig':self.test_sig,
            'mutant_rate': self.mutant_rate,
            'mut_rate_prob': self.mutant_rate_prob,
            'mut_rate_likelihood':self.mutant_rate_likelihood,
            'mut_rate_vaf':self.mutant_rate_vaf,
            'outlier_clusters': self.outlier_clusters,
            'outlier_vaf': self.outlier_vaf,
            'outlier_MI_p': handle_p_value_log10(self.outlier_moranI_pvalue),
            'outlier_MI_s': self.outlier_moranI_stat,
            'num_outlier_cluster': self.num_outlier_cluster,

            "ref_ins_major_prop": self.ref_ins_major_prop,
            # 'ref_querypos_list':self.ref_querypos_list,
            'ref_querypos_num':self.ref_querypos_num,
            # 'alt_querypos_list':self.alt_querypos_list,
            'alt_querypos_num':self.alt_querypos_num,

            'left_ref_querypos_num_remove_clip':self.left_ref_querypos_num_remove_clip,
            'left_alt_querypos_num_remove_clip':self.left_alt_querypos_num_remove_clip,

            'right_ref_querypos_num_remove_clip':self.right_ref_querypos_num_remove_clip,
            'right_alt_querypos_num_remove_clip':self.right_alt_querypos_num_remove_clip,

            'fref':self.fref,
            'falt':self.falt,

            'muts_in_cluster_p': handle_p_value_log10(self.muts_in_cluster_p),
            'muts_in_cluster_s': self.muts_in_cluster_s
        }

        for key in back_dict.keys():
            value_to_str=back_dict[key]
            if type(back_dict[key])==list:
                value_to_str=",".join([str(i) for i in back_dict[key]])
                if value_to_str=="":
                    value_to_str="no"
            if back_dict[key]=="":
                value_to_str="no"

            back_dict[key]=value_to_str

        return back_dict
    
    def test_values(self):
        test_dict={
        'identifier': self.identifier,
        # 'alt_SpotNum': self.alt_SpotNum,
        # 'dp': self.dp,
        # "ref_ins_num": self.ref_ins_num,
        # "alt_ins_num": self.alt_ins_num,
        # "ref_ins_length": self.ref_ins_length,
        # "alt_ins_length": self.alt_ins_length,
        # "ref_ins_distance": self.ref_ins_distance,
        # "alt_ins_distance": self.alt_ins_distance,

        # "ref_del_num": self.ref_del_num,
        # "alt_del_num": self.alt_del_num,
        # "ref_del_length": self.ref_del_length,
        # "alt_del_length": self.alt_del_length,
        # "ref_del_distance": self.ref_del_distance,
        # "alt_del_distance": self.alt_del_distance,

        # "ref_ins_prop": self.ref_ins_prop,
        # "ref_del_prop": self.ref_del_prop,
        # "ref_ins_major_prop": self.ref_ins_major_prop,
        # "ref_del_major_prop": self.ref_del_major_prop,
        # "alt_ins_prop":self.alt_ins_prop,
        # "alt_del_prop": self.alt_del_prop,
        # "alt_ins_major_prop":self.alt_ins_major_prop,
        # "alt_del_major_prop":self.alt_del_major_prop, 
        'combine_nearest_phase_haplotype':self.combine_nearest_phase_haplotype,
        'combine_nearest_info_mutant_prop':self.combine_nearest_info_mutant_prop,
        'combine_nearest_discordant_prop':self.combine_nearest_discordant_prop,
        'combine_nearest_phase_distance':self.combine_nearest_phase_distance,

        'no_combine_most_phase_haplotype':self.no_combine_most_phase_haplotype,
        'no_combine_most_info_mutant_prop':self.no_combine_most_info_mutant_prop,
        'no_combine_most_discordant_prop':self.no_combine_most_discordant_prop,
        'no_combine_most_phase_distance':self.no_combine_most_phase_distance,

        'no_combine_nearest_phase_haplotype':self.no_combine_nearest_phase_haplotype,
        'no_combine_nearest_info_mutant_prop':self.no_combine_nearest_info_mutant_prop,
        'no_combine_nearest_discordant_prop':self.no_combine_nearest_discordant_prop,
        'no_combine_nearest_phase_distance':self.no_combine_nearest_phase_distance,

        'combine_most_phase_haplotype':self.combine_most_phase_haplotype,
        'combine_most_info_mutant_prop':self.combine_most_info_mutant_prop,
        'combine_most_discordant_prop':self.combine_most_discordant_prop,
        'combine_most_phase_distance':self.combine_most_phase_distance

        # 'dp': self.dp,
        # 'alt_dp':self.alt_dp,
        # "reverse_dp": self.reverse_dp,
        # "forward_dp": self.forward_dp,
        # "alt_reverse_dp": self.alt_reverse_dp,
        # "alt_forward_dp": self.alt_forward_dp,
        # "major_read_strand": self.major_read_strand,
        # 'mean_distance_to_end':self.mean_distance_to_end,
        # 'RNAMutationType':self.RNAMutationType,
        # 'DNAMutationType':self.DNAMutationType,
        # 'alt_UMI_consistence_prop': self.alt_UMI_consistence_prop,
        # 'alt_consistence_hard_prop': self.alt_consistence_hard_prop,
        # 'alt_consistence_soft_prop': self.alt_consistence_soft_prop,
        # 'read_number_p':self.read_number_p,
        # 'ref_read_number_perUMI_median':self.ref_read_number_perUMI_median,
        # 'alt_read_number_perUMI_median':self.alt_read_number_perUMI_median,
        # 'ref_read_number_perUMI_max':self.ref_read_number_perUMI_max,
        # 'alt_read_number_perUMI_max':self.alt_read_number_perUMI_max,
        # 'vaf_cluster_mean':self.vaf_cluster_mean,
        # 'vaf_cluster_std_dev':self.vaf_cluster_std_dev,
        # spatial related
        # 'r2': self.r2,
        # 'wilcoxon_p':handle_p_value_log10(self.wilcoxon_p),
        # 'KS_p': handle_p_value_log10(self.KS_p),
        # 'MI_p': handle_p_value_log10(self.MI_p),
        # 'sf_test_sig':self.test_sig,
        # 'mut_rate_prob': self.mutant_rate_prob,
        # 'mut_rate_likelihood':self.mutant_rate_likelihood,
        # 'mut_rate_vaf':self.mutant_rate_vaf,
        # 'outlier_clusters': self.outlier_clusters,
        # 'outlier_vaf': self.outlier_vaf,
        # 'outlier_MI': self.outlier_moranI_pvalue,
        # 'num_outlier_cluster': self.num_outlier_cluster
        }
        # print(test_dict)
        for key in test_dict.keys():
            if test_dict[key]=="" or test_dict[key]==[]:
                test_dict[key]="no"
            elif type(test_dict[key])==list:
                count_dict = Counter(test_dict[key])
                value_to_str=",".join([str(key)+":"+str(count_dict[key]) for key in count_dict.keys()])
                test_dict[key]=value_to_str
                # if set(test_dict[key])=={0}:
                #     test_dict[key]="no"
                # else:
                #     value_to_str=",".join([str(i).strip() for i in test_dict[key]])
                    # test_dict[key]=value_to_str
            elif type(test_dict[key])==str:
                values=test_dict[key].strip()
                test_dict[key]=values
    
        return test_dict


def get_script_path():

    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)
    return script_dir


def combine_phase_sf(phase_result_file,sf_file,compare_pl_path,tmp_dir,index_in_query="0 1 2 3"):
    ## using pandas will cost more memory and more time
    # df1 = pd.read_csv(phase_result_file,header=NULL)
    # df1.columns=["chrom","pos","ref","alt","alt_origin","phase_site","gene_name","total_DP","phased_DP","phased_mutDP",\
    #              "haplotype","annotation","detailed_infomation"]
    
    # df2 = pd.read_csv(sf_file)
    # df2.columns=[]
    # merged_df = pd.merge(df1, df2, on=['Col1', 'Col2', 'Col3', 'Col4'])

    #script_dir=get_script_path()
    #compare_pl_path=os.path.join(os.path.dir(script_dir),"others","compare_files.pl")
    # print(compare_pl_path)
    tmp_filename=str(uuid.uuid4())
    compare_file=os.path.join(tmp_dir,tmp_filename)
    command=f"perl {compare_pl_path} {phase_result_file} {sf_file} {index_in_query} >{compare_file}"
    result=subprocess.run(command,shell=True,check=True)

    if result.returncode!=0:
        print(f"Something wrong when run command: {command}")
        sys.exit()

    return compare_file


def handle_h5_file(spaceranger_result_dir, count_file='raw_feature_bc_matrix.h5'):
    if spaceranger_result_dir=="":
        empty_df=""
        adata=None
    else:
        adata=read_h5ad_file_2(spaceranger_result_dir,count_file=count_file)
        if adata is None:
            empty_df=""
        else:
            spots = adata.obs_names
            empty_df=pd.DataFrame(index=spots)
            empty_df['pos_x'] = adata.obs['array_row']
            empty_df['pos_y'] = adata.obs['array_col']

    return empty_df, adata


def handle_h5ad_file(h5ad_file):
    if h5ad_file=="":
        empty_df=""
        adata=None
    else:
        import scanpy as sc
        adata=sc.read_h5ad(h5ad_file)
        if adata is None:
            empty_df=""
        else:
            spots = adata.obs_names
            empty_df=pd.DataFrame(index=spots)
            empty_df['pos_x'] = adata.obs['x_array']
            empty_df['pos_y'] = adata.obs['y_array']

    return empty_df, adata


def get_intergration_from_identifier_and_file(compare_pl_path,identifier,query_file,index_in_query="0 1 2 3"):
    """
    This function is used to get the intergration result from identifier and one special query file.
    Note: The top 4 columns of query_file must be "chrom pos ref alt"
    
    """
    perl_script=compare_pl_path
    identifier_line="\t".join(identifier.strip().split("_"))

    command=f"echo -e \"{identifier_line}\" |perl {perl_script} - {query_file} {index_in_query}|sort -u"
    try:
        result=subprocess.check_output(command,text=True,shell=True)
        return result
    except:
        print(f"Something wrong when run the command: {command}")
        return ""


def get_gene_expression(bam_file, gtf_file,out_file,htseq_path=""):
    """
    output:
    ENSG00000000003.16      974
    ENSG00000000005.6       0
    ENSG00000000419.14      1775
    ENSG00000000457.14      307
    ENSG00000000460.17      19
    ENSG00000000938.13      260
    """
    #htseq-count -f bam -r pos -s no -t exon -i gene_id IN_filter.bam /storage/douyanmeiLab/yangzhirui/Reference/gencode.v44.basic.annotation.gtf > htseq_counts.txt
    gene_count_file=out_file
    htseq_count_path=os.path.join(htseq_path,"htseq-count")
    command=f"{htseq_count_path} -f bam -r pos -s no -t exon -i gene_id {bam_file} {gtf_file} > {gene_count_file}"
    result=subprocess.run(command,text=True,shell=True,check=True)
    if result.returncode!=0:
        print(f"Something wrong when run the command: {command}")

    return gene_count_file


def get_total_gene_count(gene_count_file):
    command="awk '{ sum += $2 } END { print sum }' %s" % (gene_count_file)
    if gene_count_file!="":
        try:
            result=subprocess.check_output(command,text=True,shell=True)
            result=int(result)
        except:
            print(f"Something wrong when run the command: {command}")
            result=0
    else:
        return ""

    return result


def get_info_by_bedtools(tmp_mappbablity_file,identifier):
    chrom,pos,_,_=handle_posname(identifier)
    identifier_line="\t".join([chrom,str(pos),str(pos)])
    command=f"echo -e \"{identifier_line}\" |bedtools intersect -a {tmp_mappbablity_file} -b - "

    try:
        result=subprocess.check_output(command,shell=True,text=True)
        return result
    except:
        print(f"Something wrong when run {command}!")
        return ""


def get_info_by_grep(aim_file,query_id):
    """
    support: 
    1. knownGene_file and transcript_id
    download from ucsc goldenpath(http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz)
    or to make sure the file is created from same version of gencode, it can also be downloaded from (https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1844311714_FYTa2LuA8cLSXgDCy68JSLh9vPW4&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_regionType=genome&position=chr2%3A25%2C160%2C915-25%2C168%2C903&hgta_outputType=primaryTable&hgta_outFileName=gencodeV44_knownCDS.tsv)
    Output:
    #name   chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        proteinID       alignID
    ENST00000456328.2       chr1    +       11868   14409   11868   11868   3       11868,12612,13220,      12227,12721,14409,              uc286dmu.1
    
    2. gene_count_file and gene_id
    created by htseq-count
    Output:
    ENSG00000000003.16      974
    ENSG00000000005.6       0
    """
    if aim_file.strip().split(".")[-1]=="gz":
        cat="zcat"
    else:
        cat="cat"
    
    command=f"{cat} {aim_file} |grep {query_id} - "
    try:
        result=subprocess.check_output(command,shell=True,text=True)
        return result
    except:
        return ""

def sort_gff_file(gff3_file,tmpdir):
    #awk '$3=="exon"' /storage/douyanmeiLab/yangzhirui/Reference/gencode.v44.annotation.sort.gff3|sort -k1,1 -k4,4n -

    if gff3_file.split(".")[-1]!="gz":
        exon_sort_gff_file=os.path.join(tmpdir,os.path.basename(gff3_file)+".exon.sort.gff")
        command=f"awk '$3==\"exon\"' {gff3_file}|sort -k1,1 -k4,4n - > {exon_sort_gff_file}"
    else:
        exon_sort_gff_file=os.path.join(tmpdir,os.path.basename(gff3_file).split(".gz")[0]+".exon.sort.gff")
        command=f"gzip -d {gff3_file}|awk '$3==\"exon\"' - |sort -k1,1 -k4,4n - > {exon_sort_gff_file}"

    if os.path.exists(exon_sort_gff_file):
        pass
    else:
        result=subprocess.run(command,text=True,shell=True,check=True)
        if result.returncode!=0:
            print(f"Something wrong when run the command: {command}")

    return exon_sort_gff_file


def get_mutation_bed_file(tmp_dir,mutation_file,sep="\t"):
    tmp_index=str(uuid.uuid4())
    mutation_bed_file=os.path.join(tmp_dir, tmp_index+".mutation.tmp.bed")
    result=subprocess.run("awk -F\"%s\" '{print $1, $2, $2, $3, $4}' OFS=\"\t\" %s > %s" % (sep,mutation_file, mutation_bed_file),shell=True)
    if result.returncode!=0:
        print("Something wrong!")
    
    return mutation_bed_file


def handle_mappbablity_file(tmp_dir,mappbablity_file,mutation_bed):
    """
    mutation_bed: should be the bed format
    mappbablity_file: downloaded from https://bismap.hoffmanlab.org/, recommand umap, k=24, Multi-read
    download link: https://bismap.hoffmanlab.org/raw/hg38/k24.umap.bedgraph.gz
    """
    if mappbablity_file.split(".")[-1]=="gz":
        print(f"please check the input mappbablity file {mappbablity_file}! the compressed format is not supported.")
        print(f"=== Try to decompress {mappbablity_file} ===")
        command=f"gzip -f {mappbablity_file}"
        result=subprocess.run(command,shell=True)
        if result.returncode!=0:
            print(f"Something wrong when run {command}!")
            return ""
        else:
            mappbablity_file=mappbablity_file[: len(mappbablity_file)-3]
    elif mappbablity_file=="":
        return ""

    tmp_index=str(uuid.uuid4())
    tmp_intersect_file=os.path.join(tmp_dir,tmp_index+".used_mappablity_file.tmp")
    intersect_command=f"bedtools intersect -a {mappbablity_file} -b {mutation_bed} > {tmp_intersect_file}" 
    result2=subprocess.run(intersect_command,shell=True)
    if result2.returncode!=0:
        print(f"Something wrong when run {intersect_command}!")
    
    return tmp_intersect_file


def handle_gff_info(input_string):
    key_value_pairs = input_string.split(';')
    info_dict = {}

    for pair in key_value_pairs:
        if '=' in pair:
            key, value = pair.split('=')
            info_dict[key] = value

    return info_dict


def get_context_from_reference(reference_fasta,major_read_strand,identifier,previous_base=5):
    base=dict()
    base['A']='T'
    base['T']='A'
    base['G']='C'
    base['C']='G'
    fai_file=reference_fasta+".fai"
    if os.path.exists(fai_file):
        chr_size=get_chr_size(fai_file)
    else:
        print(f"Cannot find {fai_file}! Try to index it...")
        command=f"samtools faidx {reference_fasta}"
        result=subprocess.run(command,shell=True)
        if result.returncode!=0:
            print(f"Cannot index the {reference_fasta}! pass!")
            return ""
        else:
            chr_size=get_chr_size(fai_file)
 
    reference = Fasta(reference_fasta)
    chrom,pos,ref,alt=handle_posname(identifier)
    context_20bp=str(reference[chrom][max(1,int(pos)-11):min(int(pos)+10,int(chr_size[chrom]))])
    GCcontent=(context_20bp.count('G')+context_20bp.count('C'))/len(context_20bp)
    # context=reference[chrom][int(pos)-2:int(pos)+1]
    # context2=(base[str(reference[chrom][int(pos)-2:int(pos)-1])]+base[str(reference[chrom][int(pos)-1:int(pos)])]+base[str(reference[chrom][int(pos):int(pos)+1])])[::-1]

    ref=ref.split(",")[0]
    up_base=str(reference[chrom][pos-1-1])
    down_base=str(reference[chrom][pos+1-1])
    if ref in "CT":
        sub_name=ref+">"+alt
        dna_stat_name=up_base+"["+sub_name+"]"+down_base
    else:
        sub_name=base[ref]+">"+base[alt]
        dna_stat_name=base[down_base]+"["+ sub_name + "]"+ base[up_base]
    
    if major_read_strand=="+":
        sub_name=ref+">"+alt
        rna_stat_name=up_base+"["+sub_name+"]"+down_base
    else:
        sub_name=base[ref]+">"+base[alt]
        rna_stat_name=base[down_base]+"["+ sub_name + "]"+ base[up_base]
    
    DNAMutationType=dna_stat_name
    RNAMutationType=rna_stat_name

    # is followed by ploy alt reads?
    if major_read_strand =="-":
        previous_bases=str(reference[chrom][max(1,int(pos)-1-previous_base):int(pos)-1])
    else:
        previous_bases=str(reference[chrom][int(pos):min(int(pos)+previous_base,int(chr_size[chrom]))])
    if str(alt)*previous_base==previous_bases:
        # print(str(alt)*3)
        equal_to_previous_bases=True
    else:
        equal_to_previous_bases=False

    changed_context_bases=str(reference[chrom][max(1,int(pos)-1-previous_base):int(pos)-1])+alt+str(reference[chrom][int(pos):min(int(pos)+previous_base,int(chr_size[chrom]))])
    cause_ploy_alt=str(alt)*(previous_base+1) in changed_context_bases

    del reference
    return GCcontent, DNAMutationType, RNAMutationType, equal_to_previous_bases, cause_ploy_alt


def handle_bam_file(bam_file,chrom,pos,ref,alt,CBtag,UBtag,readLen=120):
    exist_CB=[]
    result_dict={"A":defaultdict(list), "T":defaultdict(list), "C":defaultdict(list), "G":defaultdict(list), "del":defaultdict(list)}
    for geno in "ATCG":
        result_dict[geno]["dp"]=0
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
        try:
            CB=read.get_tag(CBtag).strip()
            UB=read.get_tag(UBtag).strip()
            Name=CB+"_"+UB ## In this version, the consensus info will be ignored
        except:
            # barcode_low_quality += 1
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

                else:
                    distance_to_end=left_boundary/readLen
                    result_dict[geno]["forward_dp"]+=1
                    distance_to_end_remove_clip=left_boundary_remove_clip/len(cut_pos)

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

                UMI_name = str(UB)
                barcode_name = str(CB)

                if barcode_name not in site_barcode_UMI_dict.keys():
                    site_barcode_UMI_dict[barcode_name]=defaultdict(dict)

                if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                    site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                    site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}

                site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
                # site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
    
    for barcode in site_barcode_UMI_dict.keys():
        read_have_alt=False
        read_number_per_spot=0
        # UMI_dp+=len(site_barcode_UMI_dict[barcode].keys())
        for UMI in site_barcode_UMI_dict[barcode]:             
            count_dict=site_barcode_UMI_dict[barcode][UMI]["count"]
            # quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
            # phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
            # candidate_allele,phred=get_most_candidate_allele(phred_dict,ref)
            threshold=1
            norm_count=check_UMIconsistence_for_each_geno(count_dict,threshold)
            if norm_count!=[]:
                for geno,prop in zip("ATCG",norm_count):
                    result_dict[geno]["UMI_consistence_prop"].append(prop)

            for geno in "ATCG":
                if site_barcode_UMI_dict[barcode][UMI]["count"][geno]!=0:
                    result_dict[geno]["GenoSpotNum"]+=1
                    # print([count_dict["A"],count_dict["T"],count_dict["C"],count_dict["G"]])
                    result_dict[geno]["read_number_per_UMI"].append(count_dict[geno])
                    read_number_per_spot+=site_barcode_UMI_dict[barcode][UMI]["count"][geno]

            if site_barcode_UMI_dict[barcode][UMI]["count"][alt]!=0:
                read_have_alt=True
        
        if read_have_alt==True:
            result_dict[alt]["read_number_per_spot"].append(read_number_per_spot)
        else:
            result_dict[ref]["read_number_per_spot"].append(read_number_per_spot)

    del in_bam_read
    return result_dict, dp


def grep_dp_from_barcode_dir1(barcode_file):
    command_ref="awk '$11==0 {print $6}' %s" % (barcode_file)
    try:
        result_ref=subprocess.check_output(command_ref,text=True,shell=True)
        ref_list=[int(float(i)) for i in result_ref.strip().split("\n")]
    except:
        ref_list=[]

    command_alt="awk '$11==1 {print $6}' %s" % (barcode_file)
    try:
        result_alt=subprocess.check_output(command_alt,text=True,shell=True)
        alt_list=[int(float(i)) for i in result_alt.strip().split("\n")]
    except:
        alt_list=[]   

    return ref_list,alt_list


def grep_dp_from_barcode_dir2(barcode_file):
    command_ref="awk '$2==0 {print $1}' %s" % (barcode_file)
    try:
        result_ref=subprocess.check_output(command_ref,text=True,shell=True)
        ref_list=[int(float(i)) for i in result_ref.strip().split("\n")]
    except:
        ref_list=[]

    command_alt="awk '$2==1 {print $1}' %s" % (barcode_file)
    try:
        result_alt=subprocess.check_output(command_alt,text=True,shell=True)
        alt_list=[int(float(i)) for i in result_alt.strip().split("\n")]
    except:
        alt_list=[]   

    return ref_list,alt_list


def check_UMIconsistence_for_each_geno(count_dict,threshold=1):
    '''
    This function is used to count the consistence or not for each geno and each dict
    Version1: we want to contaion those info: A:9,T:1. Both geno A and T will be counted as 1 UMI inconsistence
    '''
    # comsistence=defaultdict(float);not_consistence=defaultdict(float)
    # # count_nonzero = sum(1 for allele in count_dict if count_dict[allele] > 0)
    # count_nonzero_allele = [allele for allele in count_dict if count_dict[allele] >=threshold]
    # if len(count_nonzero_allele)!=1:
    #     for allele in count_nonzero_allele:
    #         not_consistence[allele]+=1
    # else:
    #     comsistence[count_nonzero_allele[0]]+=1

    UMI_DP=sum(count_dict.values())
    if UMI_DP>=threshold:
        norm_count=[count_dict[geno]/UMI_DP for geno in "ATCG"]
        return norm_count
    else:
        return []
    

def hanle_bam_dp(bam_file,chrom,pos,identifier):
    in_bam_read=pysam.AlignmentFile(bam_file, "rb")
    pos_index = pos-1
    dp_far=defaultdict(list)
    dp_near=defaultdict(list)
    #for pileupcolumn in a.pileup(str(chrom), max(0,int(pos_index)-2000), min(int(pos_index)+2000,int(chr_sizes[str(chr)])), max_depth=8000):
    for pileupcolumn in in_bam_read.pileup(str(chrom), max(0,int(pos_index)-2000), int(pos_index)+2000, max_depth=8000):
        if pileupcolumn.pos==pos-2000:
            dp_far[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos-1500:
            dp_far[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos-1000:
            dp_far[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos-500:
            dp_far[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos+500:
            dp_far[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos+1000:
            dp_far[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos+1500:
            dp_far[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos+2000:
            dp_far[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos-200:
            dp_near[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos-100:
            dp_near[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos-50:
            dp_near[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos-1:
            dp_near[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos+50:
            dp_near[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos+100:
            dp_near[identifier].append(pileupcolumn.n)
        elif pileupcolumn.pos==pos+200:
            dp_near[identifier].append(pileupcolumn.n)
        
    del in_bam_read,pos_index,pos
    return dp_far,dp_near


def handle_perline(sample,
                   CBtag,
                   UBtag,
                   tmpdir,
                   mode,
                   compare_pl_path,
                   reference_fasta,
                   grep_sf_info,
                   sf_file,
                   bam_file,
                   ind_count_file,
                   ind_geno_file,
                   gene_count_file,
                   gff3_file, 
                   knownGene_file,
                   total_gene_count, 
                   combine_phase_file,
                   no_combine_phase_file,
                   empty_df,adata,
                   used_tmp_mappbablity_file, 
                   annovar_annotaion_file, 
                   vaf_cluster_file,
                   readLen,
                   barcode_dir,
                   prior,
                   line):
    """
    This function is used to get the features for per site

    """
    if line[0]=="#":
        return {}

    if grep_sf_info==1: # input is sf info, so each line has sf info
        sline=line.strip().split("\t")
        identifier="_".join(sline[0:4])
        only_pos_identifier="\t".join([sline[0],sline[1],sline[0],sline[1]])
        mut_chrom,mut_pos,mut_ref,mut_alt=handle_posname(identifier)
        mutation_features=Features(identifier)
        mutation_features.readLen=readLen

        #Feature(spatial): KS_posteior_p MI_posteior_p KS_mutant_p MI_mutant_p MI_outlier_cluster_p num_outlier_cluster
        sf_line=sline
        mutation_features.add_info_from_sf(sf_line)

    else: # input is mutation info, so sf info show be extracted from sf file
        sline=line.strip().split("\t")
        identifier=sline[0]
        mut_chrom,mut_pos,mut_ref,mut_alt=handle_posname(identifier)
        only_pos_identifier="\t".join([mut_chrom,str(mut_pos),mut_chrom,str(mut_pos)])
        mutation_features=Features(identifier)
        mutation_features.readLen=readLen

        if sf_file!="":
            sf_line=get_intergration_from_identifier_and_file(compare_pl_path,identifier,sf_file)
            if sf_line!="":
                sf_line=sf_line.strip().split("\t")
                mutation_features.add_info_from_sf(sf_line)

    print("running for ", identifier)
    #Feature(phase):discordant distancePhase
    if combine_phase_file!="":
        combine_phase_result=get_intergration_from_identifier_and_file(compare_pl_path,identifier, combine_phase_file)
        if combine_phase_result!="":
            mutation_features.add_info_from_phase(combine_phase_result,phase_type="combine")
        # for results in combine_phase_result.strip().split("\n"):
        #     if results!="":
        #         mutation_features.add_info_from_phase(results)
    if no_combine_phase_file!="":
        no_combine_phase_result=get_intergration_from_identifier_and_file(compare_pl_path,identifier, no_combine_phase_file)
        if no_combine_phase_result!="":
            mutation_features.add_info_from_phase(no_combine_phase_result,phase_type="no_combine")
        # for results in no_combine_phase_result.strip().split("\n"):
        #     if results!="":
        #         mutation_features.add_info_from_phase(results)
    #Feature: read-level features from bam (In this version, the gene information will grep from bam file)
    if bam_file!="":
        read_info_dict,dp=handle_bam_file(bam_file,mut_chrom,mut_pos,mut_ref[0],mut_alt,CBtag,UBtag,readLen)
        # print(read_info_dict)
        if dp!=0:
            mutation_features.add_read_info(read_info_dict,dp)
        # print(mutation_features)

    # add UMI_number_per_spot
    if barcode_dir!="":
        barcode_file1=os.path.join(barcode_dir,sample+"."+identifier+".barcode.mutinfo.txt")
        barcode_file2=os.path.join(barcode_dir,identifier+".barcode.mutinfo.txt")
        
        if os.path.exists(barcode_file1):
            ref_spot_dp_list, alt_spot_dp_list=grep_dp_from_barcode_dir1(barcode_file1)
            mutation_features.add_info_from_barcode_dir(ref_spot_dp_list, alt_spot_dp_list)

        elif os.path.exists(barcode_file2):
            ref_spot_dp_list, alt_spot_dp_list=grep_dp_from_barcode_dir2(barcode_file2)
            mutation_features.add_info_from_barcode_dir(ref_spot_dp_list, alt_spot_dp_list)



    #Feature: distanceExon genetype GeneNumber
    if gff3_file!="":
        # print("add gff3_info for", identifier)
        # exon_sort_gff_file=sort_gff_file(gff3_file,tmpdir)    #bedtools closest -b -  -a  <(head -n 1 /storage/douyanmeiLab/yangzhirui/01.Data_download/06.skin/04.Analysis/04.mutations/P6_ST_vis_rep2/spatial_lineager/P6_ST_vis_rep2.after_phasing_combine.txt.sort.bed)
        exon_sort_gff_file=gff3_file
        per_line="\t".join([str(i) for i in [mut_chrom,mut_pos,mut_pos]])
        tmp_bed_file=os.path.join(tmpdir,str(mut_chrom)+str(mut_pos)+".bed")
        tmp_bed=open(tmp_bed_file,"w")
        tmp_bed.write(per_line)
        tmp_bed.close()
        run_type_gff="command"

        if run_type_gff=="pybedtool":
            mutation_Bed=BedTool(tmp_bed_file)
            gff_Bed=BedTool(exon_sort_gff_file)
            mutation_closest=mutation_Bed.closest(gff_Bed)

        elif run_type_gff=="command":
            command=f"bedtools closest -b {exon_sort_gff_file}  -a {tmp_bed_file} -wa -wb"
            try:
                # print("command is ",command)
                mutation_closest=subprocess.check_output(command,text=True,shell=True)
                mutation_closest=mutation_closest.strip().split("\n")
            except:
                print(f"Something wrong when run the command: {command}")
                mutation_closest=[""]
        else:
            pass

        for line in mutation_closest:
            if line=="":
                continue

            if run_type_gff=="pybedtool":
                sline=line.fields
            elif run_type_gff=="command":
                sline=line.strip().split("\t")
            else:
                continue
            # print(sline)
            # mutation_features.add_info_from_gff(sline)
            mutation_features.add_info_from_wgEncodeGencodeExon(sline)
        mutation_features.distanceExon=min(mutation_features.distanceExon) if mutation_features.distanceExon!=[] else "no"
        # mutation_features.distanceExonIgnoreStrand=min(mutation_features.distanceExonIgnoreStrand) if mutation_features.distanceExonIgnoreStrand!=[] else "no"
        # print(mutation_features.distanceExon, mutation_features.distanceExonIgnoreStrand)
    #Feature: expression (may be deleted)
    ## count
    # if gene_count_file!="" and mutation_features.gene_id!=[]:
    #     for gene_id in mutation_features.gene_id:
    #         gene_count_info=get_info_by_grep(gene_count_file,gene_id)
    #         mutation_features.add_gene_count_info(gene_count_info)

    # ## FPKM
    # if knownGene_file!="" and mutation_features.transcript_name!=[] and total_gene_count!="":
    #     for trans_id in mutation_features.transcript_name:
    #         knownGene_info=get_info_by_grep(knownGene_file,trans_id)
    #         mutation_features.add_info_from_knownGene(knownGene_info)
    #         mutation_features.add_expression_info(total_gene_count)

    #Feature: spot gene expression variance
    # print(mutation_features.gene_name)
    if adata is not None and mutation_features.gene_name!=[]:
        for gene_name in mutation_features.gene_name:
            try:
                gene_expression = adata[:, gene_name].X
                gene_expression_dense = gene_expression.toarray() if scipy.sparse.isspmatrix(gene_expression) else gene_expression
                empty_df['gene_expression']=gene_expression_dense
                _, moran_stat, moran_pval=spatial_moran(empty_df, col="gene_expression", alpha=0.05)
                mutation_features.MI_p_expression.append(moran_pval)
                mutation_features.MI_stat_expression.append(moran_stat)

            except:
                pass

    #Feature: mappbablity score
    if used_tmp_mappbablity_file!="":
        # mappability_info=get_intergration_from_identifier_and_file(compare_pl_path,only_pos_identifier,used_tmp_mappbablity_file,index_in_query="0 1 0 1")
        mappability_info=get_info_by_bedtools(used_tmp_mappbablity_file,identifier)
        # print(mappability_info)
        for mappability_line in mappability_info.strip().split("\n"):
            if mappability_line!="":
                mutation_features.add_mappability_info(mappability_line)

    #Feature: annotation
    # if annovar_annotaion_file_all!="":
    #     # if mutation_features.phase_condition=="no":
    #         # print("The mutation origin may lost, please run phase before running")
    #     annovar_annotaion_info=get_intergration_from_identifier_and_file(compare_pl_path,identifier,annovar_annotaion_file_all,index_in_query="2 3 5 6")

    #     # annovar_annotaion_info=get_intergration_from_identifier_and_file(compare_pl_path,identifier,annovar_annotaion_file,index_in_query="0 1 3 4")
    #     for annovar_line in annovar_annotaion_info.strip().split("\n"):
    #         # print(annovar_line)
    #         if annovar_line!="":
    #             mutation_features.add_annotation_from_annovar_all(annovar_annotaion_info)
    
    # if annovar_annotaion_file_exon!="" and mutation_features.pos_anno=="exonic":
    #     # if mutation_features.phase_condition=="no":
    #         # print("The mutation origin may lost, please run phase before running")
    #     annovar_annotaion_info=get_intergration_from_identifier_and_file(compare_pl_path,identifier,annovar_annotaion_file_exon,index_in_query="4 5 7 8")
    if annovar_annotaion_file!="":
        annovar_annotaion_info=get_info_by_grep(annovar_annotaion_file,identifier)
        if annovar_annotaion_info!="":
            mutation_features.add_annotation_from_annovar(annovar_annotaion_info)

    if ind_count_file!="":
        ind_count_info=get_intergration_from_identifier_and_file(compare_pl_path,only_pos_identifier,ind_count_file,index_in_query="0 1 0 1")
        if ind_count_info!="":
            mutation_features.add_consensus_info(ind_count_info)

    if ind_geno_file!="":
        ind_geno_info=get_intergration_from_identifier_and_file(compare_pl_path,identifier,ind_geno_file,index_in_query="0 1 3 4")
        # ind_geno_info=get_intergration_from_identifier_and_file(compare_pl_path,only_pos_identifier,ind_geno_file,index_in_query="0 1 0 1")
        
        if ind_geno_info!="":
            mutation_features.add_info_from_ind_genotype(ind_geno_info)
        else:
            print(identifier)

    if vaf_cluster_file!="":
        vaf_cluster_info=get_intergration_from_identifier_and_file(compare_pl_path,identifier,vaf_cluster_file,index_in_query="0 1 3 4")
        if vaf_cluster_info!="":
            mutation_features.add_info_from_cluster_vaf(vaf_cluster_info)

    if reference_fasta!="":
        GCcontent,DNAMutationType,RNAMutationType,equal_to_previous_bases, cause_ploy_alt=get_context_from_reference(reference_fasta,mutation_features.major_read_strand,identifier)
        mutation_features.add_context_info(GCcontent,DNAMutationType,RNAMutationType, equal_to_previous_bases, cause_ploy_alt)

    if prior!="":
        prior_info=get_intergration_from_identifier_and_file(compare_pl_path,only_pos_identifier,prior,index_in_query="0 1 0 1")
        # print(prior_info)
        if prior_info!="":
            mutation_features.add_prior(prior_info)

    # Feature: exonDesOrder
    mutation_features.expand_features()
    
    if mode=="normal":
        out_dict = mutation_features.to_dict()
    elif mode=="test":
        # return mutation_features.test_values()
        out_dict = mutation_features.test_values()

    del mutation_features
    return out_dict


def main():
    print("START:",datetime.now())
    check_dir(args.outdir)
    tmp_name=str(uuid.uuid4())
    tmpdir=os.path.join(str(args.outdir),"tmp"+tmp_name)
    check_dir(tmpdir)
    current_directory = os.path.dirname(os.path.abspath(__file__))
    compare_pl_path=os.path.join(os.path.dirname(current_directory),"others/compare_files.pl")
    # gene_count_file=args.gene_count_file
    gene_count_file=""
    gff3_file=args.gff3_file
    knownGene_file=args.knownGene_file
    combine_phase_file=args.combine_phase_file
    no_combine_phase_file=args.no_combine_phase_file
    annovar_annotaion_file=args.annovar_annotaion_file
    thread=int(args.thread)

    reference_fasta=args.fasta
    CBtag=args.CBtag
    UBtag=args.UBtag
    bam_file=args.bam
    ind_count_file=args.ind_count_file
    ind_geno_file=args.ind_geno_file
    # spaceranger_result_dir=args.spaceranger_result_dir
    h5ad_path=args.h5ad
    empty_df,adata=handle_h5ad_file(h5ad_path)
    sf_file=args.sf_file
    mutation_file=args.mutations

    # total_gene_count=get_total_gene_count(gene_count_file)
    total_gene_count=""

    if mutation_file!="":
        intersect_file=mutation_file
        lines=open(mutation_file,"r").readlines()
        sep="_"
        grep_sf_info=0

    elif sf_file!="" and mutation_file=="":
        intersect_file=sf_file
        lines=open(sf_file,"r").readlines()
        sep="\t"
        grep_sf_info=1
        
    else:
        sys.exit(f"You have to offer at least one input file! Please offer a spatial feature file or a mutation file!")

    tmp_mutation_bed=get_mutation_bed_file(tmpdir,intersect_file,sep)

    mappbablity_file=args.mappbablity_file
    used_tmp_mappbablity_file=handle_mappbablity_file(tmpdir,mappbablity_file,tmp_mutation_bed)
    
    vaf_cluster_file=args.vaf_cluster_file
    readLen=args.readLen
    mode=args.mode
    sample=args.sample
    barcode_dir=args.barcode_dir
    prior=args.prior
    if mode=="normal":   
        partial_func=partial(handle_perline,
                             sample,
                            CBtag,
                            UBtag,  
                            tmpdir,
                            mode,
                            compare_pl_path,
                            reference_fasta,
                            grep_sf_info,
                            sf_file,
                            bam_file,
                            ind_count_file,
                            ind_geno_file,
                            gene_count_file,
                            gff3_file,
                            knownGene_file,
                            total_gene_count,
                            combine_phase_file,
                            no_combine_phase_file,
                            empty_df,
                            adata,
                            used_tmp_mappbablity_file,
                            annovar_annotaion_file,
                            vaf_cluster_file,
                            readLen,
                            barcode_dir,
                            prior
                            )
    
    elif mode=="test":
        partial_func=partial(handle_perline,
                             sample,
                            CBtag,
                            UBtag,  
                            tmpdir,
                            mode,
                            compare_pl_path,
                            reference_fasta,
                            grep_sf_info,
                            sf_file, #sf_file
                            bam_file, #bam_file
                            "", #ind_count_file
                            "", #ind_geno_file
                            "",
                            "",
                            "",
                            "",
                            "", # combine_phase_file,
                            "", # no_combine_phase_file,
                            "",
                            "",
                            "",
                            "", # annovar_df
                            "", # vaf_cluster_file
                            readLen,
                            "", # barcode_dir
                            "", # prior
                            )

    # print(lines)
    with multiprocessing.Pool(thread) as pool:
        results=pool.map(partial_func, lines,chunksize=5) # type: ignore # per result is [[site_info]]
    

    # for mutation_class in results:
    #     print("===========")
    #     for attr_name in dir(mutation_class):
    #         if not callable(getattr(mutation_class, attr_name)):
    #             attr_value = getattr(mutation_class, attr_name)
    #             if attr_value not in ["no",[],""] and attr_name[0]!="_":
    #             # if attr_name=="__dict__":
    #                 print(f"{attr_name}: {attr_value}")
    outname=os.path.join(args.outdir,args.outname)
    outfile=open(outname,"w")
    header=""

    for mutation_dict in results:
        # print(mutation_dict)        
        if mutation_dict is not None and mutation_dict!={}:
            if header!="" or mutation_dict=={}:
                pass
            else:
                header="\t".join(list(mutation_dict.keys()))
                outfile.write(f'#{header}\n')
            
            per_line="\t".join([str(i) for i in list(mutation_dict.values())])
            outfile.write(f'{per_line}\n')
            # for key in mutation_dict.keys():
            #     outfile.write(f'{mutation_dict[key]}\t')
            # outfile.write(f'\n')
    del results
    outfile.close()                

    print("END:",datetime.now())
    
    shutil.rmtree(tmpdir)  


def test():
    print("yes")           

## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--sample",required=True, type=str,help="sample name")
parser.add_argument("--mutations","-m", required=False,default="",help="mutation identifier list")
parser.add_argument("--sf_file","-i", required=False,default="",help="spatial feature file")
parser.add_argument("--thread", required=False,default=2,type=int, help="thread")
parser.add_argument("--CBtag", required=False,default="CB",type=str, help="CB tag")
parser.add_argument("--UBtag", required=False,default="UB",type=str, help="UB tag")
parser.add_argument("--readLen", required=False,default=120,type=int, help="read length")
parser.add_argument("--fasta", required=False,default="", help="reference fasta file")
parser.add_argument("--combine_phase_file",required=False,default="", help="combine phase_file")
parser.add_argument("--no_combine_phase_file",required=False,default="", help="no combine phase_file")
parser.add_argument("--gff3_file",required=False,default="", help="gff3_file")
parser.add_argument("--knownGene_file",required=False,default="", help="knownGene_file")
parser.add_argument("--gene_count_file",required=False,default="", help="gene_count_file")
parser.add_argument("--vaf_cluster_file",required=False,default="", help="vaf_cluster_file")
parser.add_argument("--bam",required=False,default="", help="bam_file")
parser.add_argument("--ind_count_file",required=False,default="", help="ind_count_file")
parser.add_argument("--ind_geno_file",required=False,default="", help="ind_geno_file")
parser.add_argument("--mappbablity_file",required=False,default="", help="mappbablity_file")
parser.add_argument("--annovar_annotaion_file",required=False,default="", help="annovar_annotaion_file")
parser.add_argument("--h5ad",required=False,default="", help="The h5ad file")
parser.add_argument("--spaceranger_result_dir",required=False,default="", help="The directory storing spaceranger result, such as P2/outs")
parser.add_argument("--barcode_dir",required=False,default="", help="The directory storing barcode information (under spatial feature)")
parser.add_argument("--prior",required=False,default="", help="prior file")
parser.add_argument("--outdir",required=False, default="./",help="output directory")
parser.add_argument("--outname",required=False, default="features.txt",type=str,help="output file name")
parser.add_argument("--mode",required=False,default="normal",choices=["test","normal"], help="which kind of mode you want to use?")
args = parser.parse_args()
    
if __name__ == '__main__':
    main()  
   




