import logging
from pickle import TRUE
from pkgutil import read_code
import subprocess
import pysam
import sys
import gzip
import os
import re
from pybedtools import BedTool
from collections import defaultdict
from myutils import check_dir, get_chrom_list_from_list
import vcfpy
import pandas as pd

## bam file read
def read_bam(bam_file):
    
    f1 = pysam.AlignmentFile(bam_file,"r")
    return f1


## barcode file read (for "tissue_positions.csv" file )
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
    count = 0
    for line in f.readlines():
        count+=1
        s = line.strip().split(",")
        barcode=s[0]
        if len(s)==1:
            barcode_dict[barcode] = count
        else:
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

## bed file read (for the result of mosdepth)
def handel_bed(bed_file):
    '''
    input:
    chr1    0       14048   0
    chr1    14048   14069   1
    chr1    14069   14078   2
    
    output:
    list: [("chr1",0,14048,0),...]
    '''
    pos_list=list()
    if bed_file.split(".")[-1]=="gz":
        f = gzip.open(bed_file, 'rb')
        for line in f.readlines():
            s = line.decode().strip().split()
            chr=s[0]; pos1=int(s[1]); pos2=int(s[2]); count=int(s[3])
            pos_list.append((chr,pos1,pos2,count))
    else:
        f =open(bed_file,"r")
        for line in f.readlines():
            s = line.strip().split()
            chr=s[0]; pos1=int(s[1]); pos2=int(s[2]); count=int(s[3])
            pos_list.append((chr,pos1,pos2,count))
    #print(pos_list)
    return pos_list


## vcf file read 
def handle_vcf(vcf_file):
    '''
    input:
    chrX    119811135       .       T       C,<*>   0 
    chrX    119811136       .       T       C,<*>   0 
    chrX    119811150       .       T       C,<*>   0 
    
    output:
    list: [("chrX", "119811135", "T", "C,<*>"),...]
    each site will be restored as one turple
    '''
    mutation_identifier_list=[]
    if vcf_file.split(".")[-1]=="gz":
        f=gzip.open(vcf_file,"rb")
        for line in f.readlines():
            if line[0]!="#":
                s = line.decode().strip().split()
                chrom=s[0];pos=s[1];ref=s[3];alt=s[4]
                mutation=(str(chrom),str(pos),str(ref),str(alt))
                mutation_identifier_list.append(mutation)
    else: 
        f=open(vcf_file,"r")
        for line in f.readlines():
            if line[0]!="#":
                s = line.strip().split()
                chrom=s[0];pos=s[1];ref=s[3];alt=s[4]
                if ref not in "ATCG":
                    continue
                mutation=(str(chrom),str(pos),str(ref),str(alt))
                if mutation not in mutation_identifier_list:
                    mutation_identifier_list.append(mutation)
    return mutation_identifier_list


# def handle_spot_geno_to_get_mutation_list(spot_geno_file):
#     '''
#     input:
#     chrom	pos ID	ref	alt	spot_barcode	consensus_read_count	p_mosaic	true_mutation	max_ind_geno	max_spot_geno	max_p	G_spot_max	ind_ref	ind_alt
#     chr12	52487210	.	A	.	ATCGCTGCGTGCAGCA-1	46,0,0,0	0.7014549436733754	0	mosaic	mosaic	0.7014549436733754	0/1	A	T
    
#     output: get new allele info from the "ind_ref" and "ind_alt"
#     list: [("chrX", "119811135", "T", "C,<*>"),...]
#     '''
#     mutation_identifier_list=[]
#     f=open(spot_geno_file,"r")
#     for line in f.readlines():
#         if line[0]!="#" and line[0:5]!="chrom":
#             s = line.strip().split()
#             chrom=s[0];pos=s[1];ref=s[13];alt=s[14]
#             mutation=(str(chrom),str(pos),str(ref),str(alt))
#             if mutation not in mutation_identifier_list:
#                 mutation_identifier_list.append(mutation)

#     return mutation_identifier_list


def read_file_with_multi_alt_to_identifier(file):
    '''
    lysis errors:
    chr10   5498745 5498745 T       A,G     -
    chr10   5498768 5498768 G       A,T     -
    '''
    mutation_identifier_list=[]
    f=open(file,"r")
    for line in f.readlines():
        if line[0]!="#" and line[0:5]!="chrom":
            s = line.strip().split()
            chrom=s[0];pos=s[1];ref=s[3];alts=s[4]
            for alt in alts.split(","):
                mutation="_".join([str(chrom),str(pos),str(ref),str(alt)])
                if mutation not in mutation_identifier_list:
                    mutation_identifier_list.append(mutation)

    return mutation_identifier_list


def handle_posname(pos_name):
    sitem=pos_name.split("_")
    chrom=str(sitem[0])
    pos=int(sitem[1])
    ref=str(sitem[2])
    alt=str(sitem[3])

    return chrom, pos, ref, alt

def read_mutation_list(file):
    '''
    chr12_52487210_A_T
    chr12_52487211_A_T
    '''
    mutation_identifier_list=[]
    f=open(file,"r")
    for line in f.readlines():
        if line[0]!="#" and line[0:5]!="chrom":
            mutation = line.strip()
            chrom,pos,ref,alt=handle_posname(mutation)
            for item in alt.split(","):
                mutation="_".join([str(i) for i in [chrom,pos,ref,item]])
                if mutation not in mutation_identifier_list:
                    mutation_identifier_list.append(mutation)

    return mutation_identifier_list


def handle_spot_geno_to_get_mutation_list(spot_geno_file):
    '''
    input:    
    #chrom  pos     ID      germline        mutant  cluster spot_barcode    consensus_read_count    l_germline      l_mosaic        max_spot_geno   G_spot_max      depth   vaf     p_mosaic
    chr12   52487210        .       A       T       1       ATCGCTGCGTGCAGCA-1      46,0,0,0        0.7874611700926878      0.21253882990731235     germline        0/0     46      0.0     0.3685326340225939
    
    output: get new allele info from the "ind_ref" and "ind_alt"
    list: [("chrX", "119811135", "T", "C,<*>"),...]
    '''
    mutation_identifier_list=[]
    f=open(spot_geno_file,"r")
    for line in f.readlines():
        if line[0]!="#" and line[0:5]!="chrom":
            s = line.strip().split()
            chrom=s[0];pos=s[1];ref=s[3];alt=s[4]
            mutation="_".join([str(chrom),str(pos),str(ref),str(alt)])
            if mutation not in mutation_identifier_list:
                mutation_identifier_list.append(mutation)

    return mutation_identifier_list


def handel_and_split_vcf_for_bed(vcf_file,split_vcf_dir):
    '''
    input:
    chrX    119811135       .       T       C,<*>   0 
    chrX    119811136       .       T       C,<*>   0 
    chrX    119811150       .       T       C,<*>   0 
    
    output:
    a bed file stored in tmp dir
    chrX    119811134   119811135   T
    chrX    119811135   1119811136       T
    chrX    119811149   119811150       T
    '''
    out_name=os.path.basename(vcf_file).split(".")[0]
    out_file=os.path.join(split_vcf_dir,out_name+".tmp.bed")
    result1=subprocess.run("grep -v \"#\" %s|awk '{print $1,$2-1,$2,$4}' OFS='\\t'>%s" \
                            "" % (vcf_file,out_file),shell=True)
    if result1.returncode != 0:
        print(result1.stderr)
        logging.error("There are something wrong when extracting bed file from vcf file, command is:\n" \
                 "grep -v \"#\" %s|awk '{print $1,$2-1,$2,$4}' OFS='\\t'>%s" \
                 "" % (vcf_file,out_file))
    
    result2=subprocess.run("awk '{print > (\"%s/split_%s_\" $1 \".bed\")}' OFS='\\t' %s" % (split_vcf_dir,out_name,out_file),shell=True)
    if result2.returncode != 0:
        print(result2.stderr)
        logging.error(f'Please check the tmp file: {out_file}')
    
    out_file_list=[]
    for files in os.listdir(split_vcf_dir):
        # print(os.path.basename(files)[0:5])
        if os.path.basename(files)[0:5]=="split":
            out_file_list.append(os.path.join(split_vcf_dir,files))
    
    # print(out_file_list)
    out_file_dcit=get_chrom_list_from_list(out_file_list)

    return out_file_dcit



## pos file read (the simplified result from vcf file, can be replaced by vcf file)
def handle_pos_bed(bed_file):
    '''
    input:
    chr1    14623   C       A,<*>
    chr1    14653   C       T,<*>
    
    output:
    list: [("chrX", "119811135", "T", "C,<*>"),...]
    '''
    mutation_identifier_list=[]
    f=open(bed_file,"r")
    for line in f.readlines():
        s=line.strip().split()
        chrom=s[0];pos=s[1];ref=s[2];alt=s[3]
        mutation=(str(chrom),str(pos),str(ref),str(alt))
        mutation_identifier_list.append(mutation)
    
    return mutation_identifier_list


def handle_prior_bed(bed_file):
    '''
    input:
    chr1    14622   14623       A
    chr1    14652   14653       T
    
    output:
    list: ["chrX_119811135_T_N",...]
    '''
    mutation_identifier_list=[]
    f=open(bed_file,"r")
    for line in f.readlines():
        s=line.strip().split()
        chrom=s[0];pos=s[2];ref=s[3];alt="N"
        mutation="_".join([str(chrom),str(pos),str(ref),str(alt)])
        mutation_identifier_list.append(mutation)
    
    return mutation_identifier_list


def handel_total_gnomad_file(whole_gnomAD):
    '''
    # whole_gnomAD="hg38_gnomad312_genome.txt"
    # this is used to split chromosome from gnomAD file downloaded from annovar
    '''
    result=subprocess.run("awk '{print > (\"hg38_gnomad312_genome_\" $1 \".txt\")}' OFS=\"\t\" {gnomAD_file};date".format(gnomAD_file=whole_gnomAD), shell=True)
    if result==0:
        return TRUE
    else:
        print(f'Please check the gnomAD file: {whole_gnomAD}')
        sys.exit()


def read_gnomAD_file_from_annovar(chr_gnomad_file):
    chr_gnomAD=BedTool(chr_gnomad_file)
    return chr_gnomAD


def get_freq_for_mutation_complete_from_annovar(chr_gnomAD, chr, pos, ref, alt):
    '''
    # Example lines:
    10140-10150(!10145!): ACCCT!A!ACCCC
    1   chr1	10144	10145	CCCTAA	C
    2   chr1	10145	10145	A	AC	0
    3   chr1	10145	10145	A	C	5.898e-05
    ## There are several conditions
    1   deltion of A, influence freqA
    2   A to AC, an insertion of C, didn't influent freqA
    3   substitutaion of A, influence freqA
    ### if input is ref=A, alt=C (substitution):
        p_ref = 1 - p_alt - p_deletion -p_other_substitution
    ### if input is ef=A, alt=AC (insertion): 
        p_ref = 1 - p_alt
    ### if input is ref=TAACCCC, alt= T (deletion):
        p_ref = 1- p_alt
    '''
    pos_bed="\t".join([chr,pos,ref,alt])
    pos_bed_obj=BedTool(pos_bed,from_string=True)
    pos_intersect = chr_gnomAD.intersect(pos_bed_obj,u=True)
    p_other_deletion=0.0
    p_other_sub=0.0
    p_alt=0.0
    for line in pos_intersect:
        sline=line.strip().split()
        if "-" not in [sline[3],sline[4]] and int(pos) in range(int(sline[1],int(sline[2]+1))):
            # in annovar file, the pos is recorded by "chr true_pos true_pos", 
            # while in our input file, the pos is recorded by "chr true_pos",
            # so the range should be "chr true_pos true_pos+1"
            if sline[0]==chr and int(sline[1])==pos and sline[2]==ref and sline[3]==alt:
                p_alt=float(sline[5])
            elif int(sline[1])==pos-1 and len(sline[3])<len(sline[4]):
                # this may be a deletion (-alt)
                p_other_deletion += float(sline[5])
            elif int(sline[1])==pos and int(sline[2])==pos and len(sline[3])>len(sline[4]):
                # this may be an insertion (-alt)
                pass
            elif int(sline[1])==pos and int(sline[2])==pos and len(sline[3])==len(sline[4]):
                # this may be a subsitution
                p_other_sub+=float(sline[5])
            else:
                pass
        
    if len(ref)<len(alt):
       # insertion
        p_ref=1-p_alt
    else:
        p_ref = 1- p_other_deletion - p_other_sub
    
    return p_ref, p_alt


def get_freq_for_mutation_simple_from_annovar(gnomad, identifier):
    '''
    identifier:"chr1_1000_A_alt"
    '''
    identifier_list=identifier.split("_")
    #print(identifier_list)
    chr=identifier_list[0]
    pos=int(identifier_list[1])
    ref=identifier_list[2]
    if chr=="chrM":
        out_list=[]
    else:
        # print(gnomad_dict)
        if type(gnomad) is dict:
            gnomad_file=choose_chr_file(gnomad,chr)
        elif type(gnomad) is str:
            gnomad_file=gnomad

        if gnomad=={} or gnomad_file=="":
            out_list=[]
        else:
            chr_gnomAD=read_gnomAD_file_from_annovar(gnomad_file)
            ref_list=[]
            pos_bed="\t".join([chr,str(pos-1),str(pos),ref])
            pos_bed_obj=BedTool(pos_bed,from_string=True)
            pos_intersect = chr_gnomAD.intersect(pos_bed_obj,u=True) # type: ignore
            p_dict={"A":0.0,"T":0.0,"C":0.0,"G":0.0}
            
            if pos_intersect.count()==0:
                p_dict[ref]=1.0
                # print("count==0")
            else:
                # for record in pos_intersect.fetch(str(chr), int(pos)-1,int(pos)):
                #     if record.is_snv() and record.POS==int(pos):
                #         # remove other mutation type
                #         mutation=record.ALT[0].value.capitalize()
                #         p_dict[mutation]=float(record.INFO["AF"][0])
                #         ref_list.append(record.REF.capitalize())
                for line in pos_intersect:
                    sline=line.fields
                    # print(sline)
                    if len(sline[3])==1 and len(sline[4])==1 and sline[3]!="-" and sline[4]!="-" and str(sline[2])==str(pos):
                        # remove del and ins, and remove other not equal to pos
                        # print(sline)
                        try:
                            p_dict[sline[4].capitalize()]+=float(sline[5])
                        except:
                            print(sline)
                            # ['chr19', '38773083', '38773083', 'T', 'C', '.', '0.5458', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
                        ref_list.append(sline[3])

                if len(set(ref_list))==1:
                    char_ref=list(set(ref_list))[0]
                    p_ref=1 - sum([p_dict[char] for char in p_dict.keys()])
                    p_dict[char_ref]=p_ref
                    if ref !="" and ref!=char_ref:
                        # print(identifier,ref,char_ref)
                        logging.error("Maybe you used the different reference version for gnomAD file and SNP file, please check it!")

                elif len(set(ref_list))==0:
                    p_dict[ref]=1.0

                else:
                    print(set(ref_list))
                    logging.error("There may be some wrong in the input gnomAD file in {pos}".format(pos=str(chr)+"_"+str(pos)))
            
            out_list=identifier_list[0:-1]+[p_dict["A"],p_dict["T"],p_dict["C"],p_dict["G"]]
            print(out_list)
    return out_list


def read_vcf_file_by_pyvcf(vcf_file):
    reader = vcfpy.Reader.from_path(vcf_file)
    return reader


def choose_chr_file(chr_file_dict,chr):
    if chr in chr_file_dict.keys():
        gnomad_file=chr_file_dict[chr]
    else:
        if "none" in chr_file_dict.keys():
            gnomad_file=chr_file_dict["none"]
        else:
            gnomad_file=""

    return gnomad_file


def get_freq_for_mutation_simple_from_gnomad(gnomad_file, identifier):
    '''
    identifier:"chr1_1000_A_alt"
    '''
    identifier_list=identifier.split("_")
    chr=identifier_list[0]
    pos=identifier_list[1]
    ref=identifier_list[2]
    if chr=="chrM":
        out_list=[]
    else:
        if type(gnomad_file) is dict:
            gnomad_file=choose_chr_file(gnomad_file,chr)

        elif type(gnomad_file) is str:
            gnomad_file=gnomad_file

        chr_reader=read_vcf_file_by_pyvcf(gnomad_file)
        ref_list=[]
        p_dict={"A":0.0,"T":0.0,"C":0.0,"G":0.0}
        for record in chr_reader.fetch(str(chr), int(pos)-1,int(pos)):
            if record.is_snv() and record.POS==int(pos):
                # remove other mutation type
                mutation=record.ALT[0].value.capitalize()
                p_dict[mutation]=float(record.INFO["AF"][0])
                ref_list.append(record.REF.capitalize())
        
        if len(set(ref_list))==1:
            char_ref=list(set(ref_list))[0]
            p_ref=1 - sum([p_dict[char] for char in p_dict.keys()])
            p_dict[char_ref]=p_ref
            if ref !="" and ref!=char_ref:
                logging.error("Maybe you used the different reference version for gnomAD file and SNP file, please check it!")

        elif len(set(ref_list))==0:
            p_dict[ref]=1.0

        else:
            logging.error("There may be some wrong in the input gnomAD file in {pos}".format(pos=str(chr)+"_"+str(pos)))

        out_list=identifier_list[0:-1]+[p_dict["A"],p_dict["T"],p_dict["C"],p_dict["G"]]

    return out_list


def read_prior_file(file):
    # the newest bedtools version cannot read the float file
    # bed_file=BedTool(file)
    if file=="":
        df=pd.DataFrame()
    else:
        df = pd.read_csv(file,sep='\t',header=None)
        colnames=["chr","pos","ref","fA","fT","fC","fG"]
        df.columns=colnames
        df['identifier'] = df['chr']+"_"+df['pos'].astype(str)
        df=df.set_index('identifier')
    return df


## if the file is too large, use bedtools instead
def read_count_file(file):
    mutation_list=[]
    f=open(file,"r")
    for line in f.readlines():
        s=line.strip().split()
        mutation_list.append(s)
    return mutation_list


def read_count_file_df(file):
    """read count file as a dataframe"""
    colnames=["chr","pos","ID","ref","alt","spot_barcode","consensus_read_count","qA","qT","qC","qG"]
    df = pd.read_csv(file,sep='\t',header=None, names=colnames, comment = "#",keep_default_na=False)
    df['pos'] = df['pos'].astype(str)
    return df


def read_ind_posterior_file(file, index=True):
    # df = pd.read_csv(file,sep='\t',header=None)
    colnames=["chr","pos","ID","germline","mutant","cluster","spot_number","consensus_read_count",\
              "genotype","p_mosaic","Gi","vaf"]
    df = pd.read_csv(file,sep='\t',header=None, names=colnames, comment = "#",keep_default_na=False)
    df['pos'] = df['pos'].astype(str)
    if index:
        df['identifier'] = df['chr']+"_"+df['pos'].astype(str)
        df=df.set_index('identifier')
    return df


# def read_cluster_posterior_file(file, index=True):
#     colnames=["chr","pos","ID","ref","alt","cluster","spot_number","consensus_read_count","p_refhom","p_althom","p_het","p_somatic","genotype","posterior","Gi",\
#               "germline","mutant","vaf"]
#     df = pd.read_csv(file,sep='\t',header=None, names=colnames, comment = "#",keep_default_na=False)
#     df['pos'] = df['pos'].astype(str)
#     df['cluster'] = df['cluster'].astype(str)
#     if index:
#         df['identifier'] = df['chr']+"_"+df['pos']+"_"+df['cluster'].astype(str)
#         df=df.set_index('identifier')
#     return df


def read_cluster_vaf_file(file, index=True):
    #chrom	site	ID	germline mutant	cluster spot_number	consensus_read_count  vaf
    colnames=["chr","pos","ID","germline","mutant","cluster","spot_number","consensus_read_count","vaf"]
    df = pd.read_csv(file,sep='\t',header=None, names=colnames, comment = "#",keep_default_na=False)
    df['pos'] = df['pos'].astype(str)
    df['cluster'] = df['cluster'].astype(str)
    if index:
        df['identifier'] = df['chr']+"_"+df['pos']+"_"+df['cluster'].astype(str)
        df=df.set_index('identifier')
    return df


def read_cell_number_file(file):
    colnames=["spot_barcode", "cluster", "nUMI", "nREAD", "cell_num"]
    df = pd.read_csv(file,sep='\t',header=0, names=colnames, comment = "#",keep_default_na=False)
    df['cell_num'] = df['cell_num'].astype(float)
    return df


def read_spot_posterior_file(file):
    # df = pd.read_csv(file,sep='\t',header=None)
    colnames=["chr","pos","ID","germline","mutant","cluster","spot_barcode","consensus_read_count", \
              "l_germline", "l_mosaic", "max_spot_geno","G_spot_max","depth","vaf","p_mosaic"]
    df = pd.read_csv(file, sep='\t', header=None, names=colnames, comment = "#")
    # df.columns=colnames
    df['pos'] = df['pos'].astype(str)
    df['identifier'] = df['chr']+"_"+df['pos']+"_"+df['germline']+"_"+df['mutant'].astype(str)
    df = df.set_index('identifier')
    # spot_geno_file=BedTool(file)
    return df


# def read_spot_cluster_posterior_file(file):
#     colnames=["chr","pos","ID","ref","alt","cluster","spot_barcode","consensus_read_count", \
#               "l_refhom", "l_althom", "l_het", "l_somatic", \
#               "max_cluster_geno", "max_spot_geno","max_p","G_spot_max","germline","mutant","depth","vaf","p_mosaic"]
#     df = pd.read_csv(file, sep='\t', header=None, names=colnames, comment = "#")
#     # df.columns=colnames
#     df['pos'] = df['pos'].astype(str)
#     df['cluster'] = df['cluster'].astype(str)
#     df['identifier'] = df['chr']+"_"+df['pos']+"_"+df['germline']+"_"+df['mutant'].astype(str)
#     # df=df.set_index('identifier')
#     # spot_geno_file=BedTool(file)
#     return df


def read_cluster_file(file, index=True):
    df = pd.read_csv(file, sep="\t", header=None, names=['spot_barcode', 'cluster'], comment = "#",keep_default_na=False)
    df['cluster'] = df['cluster'].astype(str)
    if index:
        df=df.set_index('spot_barcode')
    return df


def read_h5ad_file_1(output_dir,data_dir,sample_name,h5_file_name,spatial_file_name,image_file_name):
    """
    method 1: read file by scanpy.10x_h5
    data_dir: store all result from cellranger/spaceranger, such as: P6_ST_vis_rep2/outs
    """
    import scanpy as sc
    import cv2


    res_dir = os.path.join(output_dir, "result")
    check_dir(res_dir)
    
    # Read in gene expression and spatial location
    h5ad_file = os.path.join(output_dir, sample_name+"_data.h5ad")
    ## Read original 10x_h5 data and save it to h5ad
    h5_file = os.path.join(data_dir, h5_file_name)
    adata = sc.read_10x_h5(h5_file)
    spatial_file = os.path.join(data_dir, spatial_file_name)
    spatial = pd.read_csv(spatial_file,sep=",",header='infer',na_filter=False,index_col=0) 
    adata.obs["in_tissue"] = spatial['in_tissue']
    adata.obs["x_array"] = spatial['array_row']
    adata.obs["y_array"] = spatial['array_col']
    adata.obs["x_pixel"] = spatial['pxl_row_in_fullres']
    adata.obs["y_pixel"] = spatial['pxl_col_in_fullres']
    ## Select captured samples
    adata = adata[adata.obs["in_tissue"]==1]
    adata.var_names = [i.upper() for i in list(adata.var_names)]
    adata.var["genename"] = adata.var.index.astype("str")
    adata.write_h5ad(h5ad_file)
        
    ## Read in hitology image
    img_file = os.path.join(data_dir, image_file_name)
    img = cv2.imread(img_file)

    # # Set coordinates
    # x_array = adata.obs["x_array"].tolist()
    # y_array = adata.obs["y_array"].tolist()
    # x_pixel = adata.obs["x_pixel"].tolist()
    # y_pixel = adata.obs["y_pixel"].tolist()

    # # Test coordinates on the image
    # img_new = img.copy()
    # for i in range(len(x_pixel)):
    #     x=x_pixel[i]
    #     y=y_pixel[i]
    #     img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0
    # map_file = os.path.join(output_dir, sample_name+"_map.jpg")
    # cv2.imwrite(map_file, img_new)

    return adata


def read_h5ad_file_2(data_dir,count_file='raw_feature_bc_matrix.h5'):
    """
    method 2: read file by scanpy.read_visium
    data_dir: store all result from cellranger/spaceranger, such as: P6_ST_vis_rep2/outs
    """

    import scanpy as sc
    import scipy
    import warnings
    
    raw_position_file=os.path.join(data_dir,"spatial/tissue_positions.csv")
    aim_position_file=os.path.join(data_dir,"spatial/tissue_positions_list.csv")
    if not os.path.exists(aim_position_file):
        command=f"ln -s {raw_position_file} {aim_position_file}"
        result=subprocess.run(command,shell=True)
        if result.returncode!=0:
            print(f"Something wrong when run command: {command}")
    
    try:
        warnings.filterwarnings('ignore', category=UserWarning)
        adata=sc.read_visium(data_dir, genome=None, count_file=count_file, library_id=None, load_images=True, source_image_path=None)
        warnings.filterwarnings('default', category=UserWarning)
    except:
        print(f"Something wrong when run command: sc.read_visium({data_dir}, count_file={count_file})")
        adata=None

    # adata=sc.read_visium(data_dir, genome=None, count_file=count_file, library_id=None, load_images=True, source_image_path=None)    
    return adata