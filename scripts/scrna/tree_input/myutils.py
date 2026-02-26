import collections
import multiprocessing
import subprocess
from pybedtools import BedTool
import os
from pickle import TRUE
import re
import sys
from math import exp, floor, log, log10,ceil,inf
import argparse
import numpy as np
from scipy import stats
from collections import defaultdict
from functools import reduce

# from module.UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele, handle_pos, handle_quality_matrix, handle_seq


def phred_2_q(phred):
    try:
        phred=int(phred)
        q=1 - (10 ** -(phred/10))
    except ValueError as e:
        q=0
        print(f"wrong phred informat: {phred}")
        print("Error:", e)
    return q


def q_2_phred(q):
    try:
        q=float(q)
        if q==1.0:
            phred=40
        else:
            phred=ceil(-log10(1-q) * 10)
    except ValueError as e:
        print(f"wrong q informat: {q}")
        print("Error:", e)
        # sys.exit()
        phred=0
    return phred



def trans(qual_list):
    if qual_list ==[]:
        qual_str="NA"
    else:
        qual_str=",".join(qual_list)
    return qual_str


def get_chrom_list_from_file(file):
    '''
    gnomad file list or one gnomad file
    gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz
    gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz
    ...
    gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz
    gnomad.genomes.v3.1.2.sites.chrY.vcf.bgz
    gnomad.genomes.v3.1.sites.chrM.vcf.bgz
    '''
    gnomad_file_list=[];chr_list=[]
    suffix=str(file).split(".")[-1]
    if suffix in ["gz","bgz","vcf"]:
        # this may be a gnomad file, not a list file
        gnomad_file_list=[file]
        chr_list=["none"]
    else:
        # this may be a gnomad list file
        with open(file ,"r")as f1:
            gnomad_file_list=[sline.strip() for sline in f1.readlines()]
        elements_all=[]
        for files in gnomad_file_list:
            elements=re.split("\.|_",str(os.path.basename(files)))
            elements_all+=elements
        chr_list=[k for k, v in dict(collections.Counter(elements_all)).items() if v == 1]

    gnomad_dict={}
    for gnomad_file, chrom in zip(gnomad_file_list, chr_list):
        gnomad_dict[str(chrom)]=gnomad_file
    return gnomad_dict


def get_chrom_list_from_list(in_list):
    '''
    the file list like:
    [file1,file2,file3]
    '''
    file_dict={}
    #print(in_list)
    if len(in_list)==1:
        chrom=re.split("\.|_",os.path.basename(in_list[0]))[0:-1][-1]
        file_dict[chrom]=in_list[0]
    else:
        elements_all=[];chr_list=[]
        for files in in_list:
            elements=re.split("\.|_",str(os.path.basename(files)))
            elements_all+=elements
        chr_list=[k for k, v in dict(collections.Counter(elements_all)).items() if v == 1]

        # print(chr_list)
        for file,chrom in zip(in_list,chr_list):
            file_dict[str(chrom)]=file
    print(file_dict)
    return file_dict


def get_process_cores(thread):
    cores_per_cpu=multiprocessing.cpu_count()
    try:
        n_cpu=int(thread)
    except:
        print(f"wrong provided cpu numbers{thread}")
        # sys.exit()
    process_number=(cores_per_cpu-1) * n_cpu

    return process_number
    

def check_dir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.makedirs(dir,exist_ok=True)
    

def check_file(file,delete=False):
    if os.path.exists(file) and delete==True:
        os.remove(file)
        return False
    elif os.path.exists(file) and delete==False:
        return True
    else:
        return False


def intersect_pos_gnomad(dict1,dict2,out_dir):
    # dict1 is pos file dict, dict2 is annovar file dict
    out_dict={}
    # print(dict1,dict2)
    for key in dict1.keys():
        if key in dict2.keys():
            out_name=os.path.join(out_dir,"short_gnomad_"+key+".vcf")
            # print(out_name)
            # print(check_file(out_name))
            if not check_file(out_name):
                chr_POS=BedTool(dict1[key])
                chr_gnomAD=BedTool(dict2[key])
                # print(key, dict1[key],dict2[key])
                chr_gnomAD.intersect(chr_POS,u=True).saveas(out_name, trackline="#track name='short gnomad file of {chrom}' color=128,0,0".format(chrom=key)) # type: ignore
            out_dict[key]=out_name
        else:
            pass
    # print(out_dict)
    return out_dict


def intersect_pos_gnomad_only_2file(file1,file2,chrom,out_dir):
    # dict1 is pos file dict, dict2 is annovar file dict

    out_name=os.path.join(out_dir,"short_gnomad_"+chrom+".vcf")
    # print(out_name)
    # print(check_file(out_name))
    if not check_file(out_name):
        chr_POS=BedTool(file1)
        chr_gnomAD=BedTool(file2)
        # print(key, dict1[key],dict2[key])
        chr_gnomAD.intersect(chr_POS,u=True).saveas(out_name, trackline="#track name='short gnomad file of {chrom}' color=128,0,0".format(chrom=chrom)) # type: ignore

    return out_name


def split_vcf(vcf_file,split_vcf_dir,work_dir,split_type,split_line=10000):
    vcf_file = os.path.abspath(vcf_file)
    if split_type=="c": # split by chromosome

        command='cd %s;cat %s |cut -f 1 |grep -v "#" |sort -u|while read C; do vcftools  --chr "${C}" --vcf %s  --stdout --recode >split.${C}.vcf; done;' \
                'cd %s'%(split_vcf_dir, vcf_file, vcf_file,work_dir)
    elif split_type=="l": # spliy by line number
        command=f"cd {split_vcf_dir};grep -v '#' {vcf_file}|split -a 4 -d --additional-suffix=.vcf -l {split_line} ; cd {work_dir}"
    else:
        print(f"wrong input for split type {split_type}!")
    
    print(command)

    result=subprocess.run(command,shell=True)
    # print(result.returncode)
    if result.returncode!=0:
        print(f"Something wrong when run the command: {command}")
    else:
        print('finish split vcf file')

    files = os.listdir(split_vcf_dir)
    prefixes = [os.path.basename(f).split(".vcf")[0] + '\n' for f in files if f.endswith(".vcf")]
    vcf_list=open(os.path.join(split_vcf_dir,"vcf_list"),"w")
    vcf_list.writelines(prefixes)


def split_bam_by_chrom(bam_file,out_dir,work_path,run_index=False):
    # check_dir(out_dir)
    # bam_name=os.path.basename(bam_file).split(".")[0]
    bam_name=os.path.basename(bam_file)
    bam_file = os.path.abspath(bam_file)
    # bam_dir=os.path.dirname(bam_name)
    outname_prefix=bam_name+".REF_"
    split_bam_list=[]; split_bam_dict={}
    # run_index=False
    print(bam_file,out_dir)
    if os.listdir(out_dir)==[]:
        run_index=True
    # else:               
    #     for file in os.listdir(out_dir):
    #         if outname_prefix not in file:
    #             print(outname_prefix,file)
    #             run_index=True
    print(run_index)
    if run_index == True:
        #result=subprocess.run("echo yes",shell=True)
        result=subprocess.run("cd %s;rm *;ln -s %s .; samtools index %s ;bamtools split -in %s -reference;cd %s" \
                              "" %  (out_dir, bam_file, bam_name, bam_name,work_path),shell=True)

        if result.returncode!=0:
            print(f'Something wrong when running code:\n bamtools split -in {bam_file} -reference')
            # sys.exit()
        else:
            print("finish split bam file")
    else:
        print("split file already exist")

    for file in os.listdir(out_dir):
        # if outname_prefix in file and file.split(".")[-1]!="bai":
        if ".REF_" in file and file.split(".")[-1]!="bai":
            split_bam_list.append(file)
        if file.split(".")[-1]!="bai" and file+".bai" not in os.listdir(out_dir):
            result=subprocess.run("cd %s; samtools index %s; cd %s" %  (out_dir, file, work_path),shell=True)
            if result.returncode!=0:
                print(f'Something wrong when running code:\n bamtools split -in {bam_file} -reference')
                # sys.exit()
    elements_all=[];chr_list=[]
    for files in split_bam_list:
        # elements=str(os.path.basename(files).replace(bam_name+".REF_","")).replace(".bam","")
        # elements_all.append(elements) 
        prefix=os.path.basename(bam_file).replace(".bam","")+".REF_"
        elements=str(os.path.basename(files).replace(bam_name+".REF_","")).replace(".bam","").replace(prefix,"") 
        elements_all.append(elements) 
    print(elements_all)
    # chr_list=[k for k, v in dict(collections.Counter(elements_all)).items() if v == 1]

    # print(chr_list)
    for file,chrom in zip(split_bam_list,elements_all):
        split_bam_dict[str(chrom)]=os.path.join(out_dir,file)
    print(split_bam_dict)
    return split_bam_dict


def str2dict(q):
    """Convert quality string to dictionary"""
    dict = {int(i.split(':')[0]):int(i.split(':')[1]) for i in q.split(',') if q !="NA"}
    return dict
    

def str2bool(v):
    """ensure boolean input"""
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    

def handle_posname(pos_name):
    sitem=pos_name.strip().split("_")
    chrom=str(sitem[0])
    pos=int(sitem[1])
    ref=str(sitem[2])
    alt=str(sitem[3])

    return chrom, pos, ref, alt

def handle_p_value_log10(p_value):
    try:
        in_type=type(p_value)
        if in_type==list:
            p_list=[]
            for p_val in p_value:
                if p_val==0 or p_val=="0":
                    p_val=1e-300
                elif p_val=="NA":
                    p_val=1
                p_list.append(log10(float(p_val)))
            return p_list
        elif in_type==str:
            if p_value=="0":
                p_value=1e-300
            elif p_value=="NA":
                p_value=1
            p=log10(float(p_value))
            return p
        
        elif in_type==float or in_type==np.float64:
            if p_value == 0.0:
                p_value=1e-300
            p=log10(float(p_value))
            return p
    except:
        print("p_val:",p_value)
        return "no"


def get_chr_size(fai_file):
    """
    /storage/douyanmeiLab/yangzhirui/Reference/Cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai
    """
    chr_sizes=dict()
    for line in open(fai_file,"r"):
        line=line.rstrip()
        fields=line.split('\t')
        chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])

    return chr_sizes


def str_to_dict(one_string):
    return_dict=collections.defaultdict(int)
    for items in one_string.strip().split(","):
        if ":" in items:
            return_dict[items.split(":")[0]]+=int(items.split(":")[1])

    return return_dict


def do_wilicox_sum_test(input_1,input_2,method="two-sided",type="dict"):
    from scipy import stats
    if type=="dict":
        input_list1=[]
        input_list2=[]
        for i in input_1.keys():
            input_list1.extend(i[0]*int(i[1]))

        for m in input_2.keys():
            input_list2.extend(m[0]*int(m[1]))    
    elif type=="list":
        input_list1=input_1
        input_list2=input_2
    else:
        print(f"wrong input in do_wilicox_sum_test: {type}")
        return "no","no"
    
    z_statistic, p_value=stats.ranksums(input_list1,input_list2,alternative=method)
    z_statistic=float(z_statistic);p_value=float(p_value)
    return z_statistic,p_value


def get_geno_from_alignment(pos,item):
    def handle_cigar(ciagr_symbol):
        '''
        ## handel cigar
        # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
        # '76M1D33M139241N11M'
        # the 1st is symbol; and the 2nd is count
        # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
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

    Name,geno,quality='','',''
    
    try:
        CB=item.get_tag("CB").strip()
        UB=item.get_tag("UB").strip()
        Name=CB+"_"+UB
        pos_index=pos-1

        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq=handle_seq(item.seq, seq_cut)
        cut_pos=handle_pos(item.get_reference_positions(), pos_cut)
        raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
        quality=item.get_forward_qualities()[raw_index]

        if pos_index in cut_pos:
            # effective_DP += 1
            geno = cut_seq[cut_pos.index(pos_index)]            
    except:
        pass

    return Name,geno,quality


def get_hard_clip_count(ciagr_symbol):
    """
    same as handle_cigar, but only used to grep hard clip
    """
    left_hardclip,right_hardclip=0,0
    for cigars, i in zip(ciagr_symbol,range(1,len(ciagr_symbol)+1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol ==5 and i==1:
            left_hardclip=count
        elif symbol==5 and i==len(ciagr_symbol):
            right_hardclip=count

    return left_hardclip,right_hardclip


def get_indel_info(indel_info_list,pos_index_in_raw_matrix):
    """
    similar as handle_pos, input is [(10,1)]. 10 is the indel position to sequence start(exclude soft/hard slip seq)
    output:
    the plus or minus of indel will display the direction between indel and mutation. The plus means indel in the left of mutation
    """
    indel_num=0;indel_length=[];indel_distance=[]
    if len(indel_info_list) == 0:
        indel_num=0
    elif len(indel_info_list) >1 and pos_index_in_raw_matrix=="in":
        indel_num=len(indel_info_list)
        indel_distance=[]
    else:
        indel_num=len(indel_info_list)
        for item in indel_info_list:
            indel_length.append(item[1])
            if pos_index_in_raw_matrix!="in":
                indel_distance.append(pos_index_in_raw_matrix-item[0])
            elif pos_index_in_raw_matrix=="in":
                indel_distance.append(0)
    # indel_length="/".join([str(i) for i in indel_length])
    # indel_distance="/".join([str(i) for i in indel_distance])
    # if indel_num!=0:
    #     print(indel_info_list,pos_index_in_raw_matrix,indel_num,indel_length,indel_distance)
    return indel_num,indel_length,indel_distance


def judge_target_pos_in_indel(pos_index,ins_or_del_info_list,get_reference_positions):
    def get_pos_list(info_list):
        for item in info_list:
            one_list=[]
            index_num=item[0];count_num=item[1]
            pos=get_reference_positions[index_num-1]
            for i in range(count_num):
                one_list.append(pos+i+1)
            if pos_index+1 in one_list:
                return one_list
            else:
                continue
        return []
    pos_list=get_pos_list(ins_or_del_info_list)
    return pos_list



def judge_pos_in_indel(ins_info_list,del_info_list,get_reference_positions):
    pos_list=[]

    def get_pos_list(info_list,pos_list):
        for item in info_list:
            index_num=item[0];count_num=item[1]
            pos=get_reference_positions[index_num-1]
            for i in range(count_num):
                pos_list.append(pos+i+1)
        return pos_list

    pos_list=get_pos_list(ins_info_list,pos_list)
    pos_list=get_pos_list(del_info_list,pos_list)

    return pos_list


# def get_deletion_info(ciagr_symbol,pos_index_in_raw_matrix):
#     del_num=0;del_length=[];del_distance=[]
#     seq_length_before=0
#     for cigars, i in zip(ciagr_symbol,range(1,len(ciagr_symbol)+1)):
#         symbol = cigars[0]
#         count = cigars[1]
#         if symbol ==0:
#             # measure  the seq length 
#             seq_length_before += count
#         elif symbol==2:
#             del_num+=1
#             del_length.append(count)
#             del_distance.append(pos_index_in_raw_matrix-seq_length_before)
#         else:
#             pass
    
#     return del_num,del_length,del_distance


def combine_info_from_cigar(cigar_symbol):
    '''
    ## combine: handel cigar;  get_hard_clip_count
    # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
    # '76M1D33M139241N11M'
    # the 1st is symbol; and the 2nd is count
    # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
    '''
    seq_length_before = 0
    pos_length_before = 0

    left_hardclip,right_hardclip=0,0

    seq_cut_start = None; seq_cut_end = None
    pos_cut=[]
    del_info=[]
    for cigars, i in zip(cigar_symbol,range(1,len(cigar_symbol)+1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol in [6,7,8]:
            # an api for handeling mapping issues "HP=X"
            print(cigar_symbol)  ## LOG
        elif symbol in [0, 1, 2, 4,5]:
            # measure the seq length 
            seq_length_before += count
            if symbol == 0:
                pos_length_before += count
            elif symbol == 4:
                # whether "S" is in this read
                if i == 1:
                    # whether the "S" is in the head or tail
                    seq_cut_start = seq_length_before
                elif i ==len(cigar_symbol):
                    seq_cut_end = seq_length_before
                else:
                    print(cigar_symbol) ## LOG
            elif symbol == 1:
                # whether the "I" is in the cigar
                pos_cut.append((pos_length_before,count))
            elif symbol == 2:
                del_info.append((pos_length_before,count))
            elif symbol==5:
                print(cigar_symbol)  ## LOG
                if i==1:
                    left_hardclip=count
                elif i==len(cigar_symbol):
                    right_hardclip=count
        else:
            pass
    seq_soft_cut = (seq_cut_start, seq_cut_end)
    seq_hard_clip = (left_hardclip,right_hardclip)

    return seq_soft_cut, pos_cut, del_info, seq_hard_clip


def list2min(value):
    """find the minimum value in the comma-separated string"""
    # check if the input is a float number
    if isinstance(value, float):
        return value
    else:
        return min(float(num) for num in value.split(','))
    
    
def list2frequent(value):
    """find the most frequent element in the comma-separated string"""
    return collections.Counter(value.split(',')).most_common(1)[0][0]


## FUNCTION: adapted from bcftools
def calc_vdb(pos_list, npos,readlen=120):
    from scipy.special import erfc
    """
    # Note well: the parameters were obtained by fitting to simulated data of
    # 100bp reads. This assumes rescaling to 100bp in bcf_call_glfgen().
    
    """
    # readlen = 100
    assert npos == readlen, "npos must be equal to readlen"

    nparam = 15
    param = np.array([
        [3, 0.079, 18], [4, 0.09, 19.8], [5, 0.1, 20.5], [6, 0.11, 21.5],
        [7, 0.125, 21.6], [8, 0.135, 22], [9, 0.14, 22.2], [10, 0.153, 22.3],
        [15, 0.19, 22.8], [20, 0.22, 23.2], [30, 0.26, 23.4], [40, 0.29, 23.5],
        [50, 0.35, 23.65], [100, 0.5, 23.7], [200, 0.7, 23.7]
    ])

    dp = 0
    mean_pos = 0
    mean_diff = 0
    for i in range(npos):
        if pos_list[i] == 0:
            continue
        dp += pos_list[i]
        mean_pos += pos_list[i] * i

    if dp < 2:
        return inf  # one or zero reads can be placed anywhere

    mean_pos /= dp
    for i in range(npos):
        if pos_list[i] == 0:
            continue
        mean_diff += pos_list[i] * abs(i - mean_pos)

    mean_diff /= dp

    ipos = int(mean_diff)  # tuned for float-to-int implicit conversion
    if dp == 2:
        return (2 * readlen - 2 * (ipos + 1) - 1) * (ipos + 1) / (readlen - 1) / (readlen * 0.5)

    if dp >= 200:
        i = nparam  # shortcut for big depths
    else:
        for i in range(nparam):
            if param[i][0] >= dp:
                break

    if i == nparam:
        # the depth is too high, go with 200x
        pscale = param[nparam - 1][1]
        pshift = param[nparam - 1][2]
    elif i > 0 and param[i][0] != dp:
        # linear interpolation of parameters
        pscale = (param[i - 1][1] + param[i][1]) * 0.5
        pshift = (param[i - 1][2] + param[i][2]) * 0.5
    else:
        pscale = param[i][1]
        pshift = param[i][2]

    return 0.5 * erfc(-(mean_diff - pshift) * pscale)


def calc_mwu_biasZ(a, b, n, left_only, do_Z):
    b_empty = all(x == 0 for x in b)
    e = l = na = nb = 0
    t = 0

    if b_empty:
        for i in reversed(range(n)):
            na += a[i]
            t += (a[i] * a[i] - 1) * a[i]  # adjustment score for ties
    else:
        for i in reversed(range(n)):
            e += a[i] * b[i]
            l += a[i] * nb  # a<b
            na += a[i]
            nb += b[i]
            p = a[i] + b[i]
            t += (p * p - 1) * p  # adjustment score for ties

    if not na or not nb:
        return float('inf')

    U = l + e * 0.5  # Mann-Whitney U score
    m = na * nb / 2.0

    # With ties adjustment
    var2 = (na * nb) / 12.0 * ((na + nb + 1) - t / ((na + nb) * (na + nb - 1)))
    if var2 <= 0:
        return 0 if do_Z else 1

    if do_Z:
        # Z-score
        return (U - m) / np.sqrt(var2)

    # U score, asymmetric for some data types
    if left_only and U > m:
        return float('inf')

    if na >= 8 or nb >= 8:
        # Normal approximation
        return np.exp(-0.5 * (U - m) ** 2 / var2)

    # Placeholder for exact calculation
    # You'll need to implement or replace mann_whitney_1947
    # return mann_whitney_1947(na, nb, U) * np.sqrt(2 * np.pi * var2)
    return "Exact calculation needed"


def logsumexp2(a, b):
    return np.log(1 + np.exp(-np.abs(a - b))) + max(a, b)


def calc_SegBias_for_one_sample(ref_dp,alt_dp):
    """
    adapted from bcftools. Here we removed the condition, which have not only one sample.
    """
    nr = alt_dp  # number of observed non-reference reads
    if not nr:
        return

    avg_dp = (ref_dp + nr)   # average depth
    M = floor(nr / avg_dp + 0.5)  # approximate number of variants sampled in the population
    M = min(max(M, 1), 1)  # Clamp M between 1 and sample number(1)
    f = M / 2.0   # allele frequency
    p = nr   # number of variant reads per sample expected if variant not real (Poisson)
    q = nr / M  # number of variant reads per sample expected if variant is real (Poisson)
    sum = 0
    log2 = log(2.0)


    oi = alt_dp  # observed number of non-ref reads
    if oi:
        tmp = logsumexp2(log(2 * (1 - f)), log(f) + oi * log2 - q)
        tmp += log(f) + oi * log(q / p) - q + p
    else:
        tmp = log(2 * f * (1 - f) * exp(-q) + f**2 * exp(-2 * q) + (1 - f)**2) + p
    sum += tmp

    return sum


def get_intergration_from_identifier_and_one_file(identifier,query_file,index_in_query="0 1 2 3"):
    """
    This function is used to get the intergration result from identifier and one special query file.
    Note: The top 4 columns of query_file must be "chrom pos ref alt"
    
    """
    current_directory = os.path.dirname(os.path.abspath(__file__))
    compare_pl_path=os.path.join(current_directory,"others/compare_files.pl")
    perl_script=compare_pl_path
    identifier_line="\t".join(identifier.strip().split("_"))

    command=f"echo -e \"{identifier_line}\" |perl {perl_script} - {query_file} {index_in_query}|sort -u"
    try:
        result=subprocess.check_output(command,text=True,shell=True)
        return result
    except:
        print(f"Something wrong when run the command: {command}")
        return ""


def closest_command(tmp_bed_file,exon_sort_gff_file):

    command=f"bedtools intersect -a {exon_sort_gff_file} -b {tmp_bed_file}"
    try:
        results=subprocess.check_output(command,text=True,shell=True)

    except:
        print(f"Something wrong when run the command: {command}")
        results=""
 
    return results


def closest_pybedtools(tmp_bed_file,exon_sort_gff_file):
    mutation_Bed=BedTool(tmp_bed_file)
    # mutation_Bed=BedTool(per_line,from_string=True)
    gff_Bed=BedTool(exon_sort_gff_file)
    mutation_closest=mutation_Bed.closest(gff_Bed)
    return mutation_closest


def round_to_nearest_bin(x,bins):
    return int(np.ceil(x / bins) * bins)


def wilcoxon_with_rbc(x, y, alternative='two-sided'):
    """ 
    Perform a Wilcoxon rank-sum test with calculating the Rank-Biserial Correlation (RBC) value.
    Parameters
    ----------
    x,y : array_like
        The data from the two samples.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis. Default is 'two-sided'.
        The following options are available:
        * 'two-sided': one of the distributions (underlying `x` or `y`) is
          stochastically greater than the other.
        * 'less': the distribution underlying `x` is stochastically less
          than the distribution underlying `y`.
        * 'greater': the distribution underlying `x` is stochastically greater
          than the distribution underlying `y`.
    Returns
    -------
    statistic : float
        The test statistic under the large-sample approximation that the
        rank sum statistic is normally distributed.
    pvalue : float
        The p-value of the test.
    rbc: float
        The Rank-Biserial Correlation (RBC) value, which is a measure of
        the strength and direction of the association between the two samples.
        It ranges from -1 to 1, where:
        - 1 indicates a perfect positive association,
        - -1 indicates a perfect negative association,
        - 0 indicates no association. 
    """
    x, y = map(np.asarray, (x, y))
    n1 = len(x)
    n2 = len(y)
    alldata = np.concatenate((x, y))
    ranked = stats.rankdata(alldata)
    x_ranks = ranked[:n1]
    R1 = np.sum(x_ranks)
    # Mann-Whitney U statistic for group x
    U = R1 - (n1 * (n1 + 1)) / 2
    # Rank Biserial Correlation (rbc)
    rbc = (2 * U) / (n1 * n2) - 1
    # statistic and p-value with large number normal approximation
    statistic, p_value = stats.ranksums(x, y, alternative=alternative)
    return statistic, p_value, rbc



def calculate_rbc_for_paired_wilcoxon(x, y):
    """ 
    Calculate the Rank-Biserial Correlation (RBC) value for the paired Wilcoxon signed-rank test.
    Parameters
    ----------
    x,y : array_like
        The data from the two samples.
    Returns
    -------
    rbc: float
        The Rank-Biserial Correlation (RBC) value, which is a measure of
        the strength and direction of the association between the two samples.
        It ranges from -1 to 1, where:
        - 1 indicates a perfect positive association,
        - -1 indicates a perfect negative association,
        - 0 indicates no association. 
    """
    from sklearn.preprocessing import StandardScaler
    # format as numpy arrays
    x, y = map(np.asarray, (x, y))
    # test whether the two arrays are of the same length
    if len(x) != len(y):
        raise ValueError("Input arrays must have the same length.")
    try:
        # standardize
        X = np.array([x, y]).T 
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)  # normalizes each column independently
        x_scaled = X_scaled[:, 0]
        y_scaled = X_scaled[:, 1] 
        # Compute RBC
        diff = x_scaled - y_scaled
        non_zero = diff != 0
        ranks = stats.rankdata(np.abs(diff[non_zero]))
        signed_ranks = ranks * np.sign(diff[non_zero])
        W_plus = signed_ranks[signed_ranks > 0].sum()
        W_minus = -signed_ranks[signed_ranks < 0].sum()
        rbc = (W_plus - W_minus) / (W_plus + W_minus)
        return rbc
    except:
        return 0
    

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
    

def handel_bam_file(bam_file,chrom,pos,ref,alt,run_type,outdir,bins,readLen=120):
    
    from collections import defaultdict
    import pysam

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

    del in_bam_read
    return result_dict,dp



## based on count and quality dict per UMI, to get all candidate allele and their phred score
def calculate_UMI_combine_phred(count_dict, quality_dict,weigh=0.5):
    all_genos=["A","T","C","G"]
    pcr_error = 1e-6
    #no_pcr_error = 1.0 - 3e-5 the reference from smcount
    no_pcr_error = (1.0 - pcr_error) ** 100 # median cycle in RNA-seqing is 100 (50-150)
    rightP = 1.0
    sumP = 0.0
    dp=sum(count_dict.values())
    proP_dict=defaultdict(lambda : 1.0)
    pcrP_dict=defaultdict(float)
    likelihood_dict=defaultdict(float)
    phred_dict=defaultdict(float)
    for geno in count_dict.keys():
        ## proP_value means no sequencing error for each geno
        # the likelihood whose allele equal to geno, here the quality is the right prob for one base
        qual_geno_list=[phred_2_q(key)**int(quality_dict[geno][key]) for key in quality_dict[geno].keys()]
        qual_geno=reduce(lambda x, y: x*y, qual_geno_list)
        proP_dict[geno]*=qual_geno
        # the likelihood whose allele not equal to geno
        for other_geno in quality_dict.keys()-set([geno]):
            other_qual_geno_list = [(1-phred_2_q(key))**int(quality_dict[other_geno][key]) for key in quality_dict[other_geno].keys()]
            if other_qual_geno_list == []:
                continue
            other_qual_geno=reduce(lambda x, y: x*y, other_qual_geno_list)
            proP_dict[geno]*=other_qual_geno
        
        ## rightP means no sequencing error, or no base calling error for all base
        rightP = rightP * qual_geno
    
    for geno in all_genos:
        ## pcrP means PCR error
        count_geno = 0 if geno not in count_dict.keys() else count_dict[geno]
        ratio = ( count_geno + 0.5) / (dp + 0.5 * 4)
        pcrP = 10.0 ** (-6.0 * ratio)
        pcrP_dict[geno]=pcrP
    
    # after obtaining [sequencing_error, no_pcr_error, no_sequencing_error, pcr_error], the likelihood of each geno will be calculate
    for geno in all_genos:
        if geno in count_dict.keys():
            base_calling_error = proP_dict[geno]
            no_base_calling_error=rightP
            pcr_error=min([pcrP_dict[char] for char in pcrP_dict.keys() if char != geno])
            likelihood_value = weigh * no_pcr_error * base_calling_error + (1-weigh) * no_base_calling_error * pcr_error 
        else:
            likelihood_value = rightP
            for char in set(all_genos) - set([geno]):
                likelihood_value *= pcrP_dict[char]
                    
        likelihood_dict[geno]=likelihood_value
        sumP += likelihood_value
    
    for geno in likelihood_dict.keys():
        phred_dict[geno] = 0 if sumP <= 0 else q_2_phred(likelihood_dict[geno] / sumP)

    return phred_dict


# following the last function, to get the most candidate allele and it's phred
def get_most_candidate_allele(phred_dict,ref_allele):
    rank_list=sorted(phred_dict.items(), key = lambda item:item[1], reverse=True)
    major_allele=rank_list[0][0]; major_allele_phred=rank_list[0][1]
    # major_allele_count=count_dict[major_allele]

    # ref_allele_count=count_dict[ref_allele]
    ref_allele_phred=phred_dict[ref_allele]

    if  major_allele != ref_allele and ref_allele_phred>=major_allele_phred:
        candidate_allele=ref_allele;phred=ref_allele_phred
    else:
        candidate_allele=major_allele;phred=major_allele_phred

    return candidate_allele,phred



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


