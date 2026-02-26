from read_file import choose_chr_file,read_bam
from myutils import q_2_phred,phred_2_q
from collections import defaultdict
from functools import reduce
# import sys
# import pandas as pd

# input: bam
# output: chr pos ID ref alt spot_num UMI_count quality_A quality_T quality_C quality_G
#           1  11  .  A   T   10      9,1,0,0     30:8,20:1   30:1    na          na


### This four function are used to get right quality and base from pysam
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


## grep the site info from bam file (the annotated line are get the statistic count information per site)
def handle_reads_per_pos_read_count(bam_handle,chrom,pos,run_type):
    '''
    input:
    bam_handel: the initial bam file handeled by pysam
    
    output:
    dict: {barcode_name: {UMI_name: {"count": 1, "quality": {"A":{30:9, 10:1}, "T":{10:1}}}}}
    '''
    # reads=bam_handle.fetch(chrom,pos,pos+1)
    reads=bam_handle.fetch(chrom,pos-1,pos,multiple_iterators=True)
    pos_index = pos-1
    # i=1
    # all_reads_count=0
    # effective_DP=0
    # mapping_reads=0;skip_reads=0
    # barcode_low_quality=0   

    site_barcode_UMI_dict={}

    for item in reads:
        # all_reads_count += 1
        try:
            # a part of reads didn't have the information of "CB", bacause the "CR" didn't pass QC
            if run_type=="visium":
                CB=item.get_tag("CB").strip()
                UB=item.get_tag("UB").strip()

                barcode_name=str(CB)
                UMI_name=str(UB)

            elif run_type=="stereo":
                Cx=str(item.get_tag("Cx"))
                Cy=str(item.get_tag("Cy"))
                UR=item.get_tag("UR").strip()

                barcode_name=Cx+"_"+Cy
                UMI_name=str(UR)

            elif run_type=="ST":
                CB=str(item.get_tag("B0"))
                UB=str(item.get_tag("B3"))

                barcode_name=str(CB)
                UMI_name=str(UB)
            else:
                # print("type",run_type)
                continue
                
        except:
            # barcode_low_quality += 1
            # print("no CB UB for ", item)
            continue
        
        try:
            item.get_reference_positions().index(pos_index)
        except:
            # skip_reads += 1
            # print("no pos_index for ", item)
            continue
        # mapping_reads += 1

        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq=handle_seq(item.seq, seq_cut)
        cut_pos=handle_pos(item.get_reference_positions(), pos_cut)

        if pos_index in cut_pos:
            # effective_DP += 1
            geno = cut_seq[cut_pos.index(pos_index)]
            if geno not in "ATCG":
                print("not ATCG ", item)
                continue

            raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
            try:
                qualities=item.get_forward_qualities()
                quality=qualities[raw_index]
            except:
                # qualities=None
                print("not qualities ", item)
                # i+=1
                continue

            if barcode_name not in site_barcode_UMI_dict.keys():
                site_barcode_UMI_dict[barcode_name]=defaultdict(dict)

            if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
 
            site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
            site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
    
    # print(i)
    #stat_info = [all_reads_count, effective_DP, barcode_low_quality, mapping_reads, skip_reads]
    #print(stat_info)
    return site_barcode_UMI_dict


def handle_reads_per_pos_read_count_and_strand(bam_handle,chrom,pos,run_type):
    '''
    input:
    bam_handel: the initial bam file handeled by pysam
    
    output:
    dict: {barcode_name: {UMI_name: {"count": 1, "quality": {"A":{30:9, 10:1}, "T":{10:1}}}}}
    '''
    # reads=bam_handle.fetch(chrom,pos,pos+1)
    reads=bam_handle.fetch(chrom,pos-1,pos,multiple_iterators=True)
    pos_index = pos-1
    # i=1
    # all_reads_count=0
    # effective_DP=0
    # mapping_reads=0;skip_reads=0
    # barcode_low_quality=0   

    site_barcode_UMI_dict={}
    reverse_dp=0
    forward_dp=0
    for item in reads:
        # all_reads_count += 1
        try:
            # a part of reads didn't have the information of "CB", bacause the "CR" didn't pass QC
            if run_type=="visium":
                CB=item.get_tag("CB").strip()
                UB=item.get_tag("UB").strip()

                barcode_name=str(CB)
                UMI_name=str(UB)

            elif run_type=="stereo":
                Cx=str(item.get_tag("Cx"))
                Cy=str(item.get_tag("Cy"))
                UR=item.get_tag("UR").strip()

                barcode_name=Cx+"_"+Cy
                UMI_name=str(UR)
            elif run_type=="ST":
                CB=str(item.get_tag("B0"))
                UB=str(item.get_tag("B3"))

                barcode_name=str(CB)
                UMI_name=str(UB)
            else:
                # print("type",run_type)
                continue
                
        except:
            # barcode_low_quality += 1
            # print("no CB UB for ", item)
            continue
        
        try:
            item.get_reference_positions().index(pos_index)
        except:
            # skip_reads += 1
            # print("no pos_index for ", item)
            continue
        # mapping_reads += 1


        if item.is_reverse in [True,"TRUE","true","True"]:
            reverse_dp+=1
        else:
            forward_dp+=1        
        

        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq=handle_seq(item.seq, seq_cut)
        cut_pos=handle_pos(item.get_reference_positions(), pos_cut)

        if pos_index in cut_pos:
            # effective_DP += 1
            geno = cut_seq[cut_pos.index(pos_index)]
            if geno not in "ATCG":
                print("not ATCG ", item)
                continue

            raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
            try:
                qualities=item.get_forward_qualities()
                quality=qualities[raw_index]
            except:
                # qualities=None
                print("not qualities ", item)
                # i+=1
                continue

            if barcode_name not in site_barcode_UMI_dict.keys():
                site_barcode_UMI_dict[barcode_name]=defaultdict(dict)

            if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
 
            site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
            site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
    
    # print(i)
    #stat_info = [all_reads_count, effective_DP, barcode_low_quality, mapping_reads, skip_reads]
    #print(stat_info)
    if reverse_dp>=forward_dp:
        major_read_strand="-"
    elif reverse_dp<forward_dp:
        major_read_strand="+"
    else:
        major_read_strand="unknown"
        
    return site_barcode_UMI_dict, major_read_strand

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


### combine UMI for individual data
def combine_UMI_bulk(site_barcode_UMI_dict,chrom,pos,ref):
    #out_df = pd.DataFrame(columns=["chr","pos","ref","alt","spot_number","consensus_read_count","read_quality_A","read_quality_T","read_quality_C","read_quality_G"])
    spot_number=len(site_barcode_UMI_dict.keys())
    all_genos=["A","T","C","G"]
    raw_count_dict=defaultdict(int)
    consensus_read_count=defaultdict(int)
    consensus_read_quality=defaultdict(dict)
    for barcode in site_barcode_UMI_dict.keys():

        for UMI in site_barcode_UMI_dict[barcode]:
            count_dict=site_barcode_UMI_dict[barcode][UMI]["count"]
            quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
            phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
            candidate_allele,phred=get_most_candidate_allele(phred_dict,ref)

            consensus_read_count[candidate_allele]+=1
            if str(phred) not in consensus_read_quality[candidate_allele].keys():
                consensus_read_quality[candidate_allele][str(phred)]=0

            consensus_read_quality[candidate_allele][str(phred)]+=1

            for geno in count_dict.keys():
                raw_count_dict[geno]+=count_dict[geno]
    # remove N allele of ref
    # if ref == "N":
    #     ref=sorted(consensus_read_count.items(), key = lambda item:item[1], reverse=True)[0][0]
    # calculate alt
    alt=",".join([alt_geno for alt_geno in consensus_read_count.keys() if alt_geno != ref])
    if alt =="":
        alt="."

    consensus_read_count_list=[]
    raw_read_count_list=[]

    consensus_read_qual_dict=defaultdict(str)
    for geno in all_genos:
        geno_count=0 if geno not in consensus_read_count.keys() else consensus_read_count[geno]
        consensus_read_count_list.append(str(geno_count))

        raw_count=0 if geno not in raw_count_dict.keys() else raw_count_dict[geno]
        raw_read_count_list.append(str(raw_count))

        if geno not in consensus_read_quality.keys():
            geno_qual_str="NA" 
        else:
            geno_qual_list=[str(phred)+":"+str(consensus_read_quality[geno][phred]) for phred in consensus_read_quality[geno].keys()]
            geno_qual_str=",".join(geno_qual_list)
        consensus_read_qual_dict[geno]=geno_qual_str

    consensus_read_count_str=",".join(consensus_read_count_list)
    raw_read_count_str=",".join(raw_read_count_list)
    site_info=[chrom, pos,".", ref, alt, spot_number, consensus_read_count_str,consensus_read_qual_dict["A"],consensus_read_qual_dict["T"],consensus_read_qual_dict["C"],consensus_read_qual_dict["G"],raw_read_count_str]
    #site_info={"chrom":chrom, "site":pos, "ID":".","ref":ref, "alt":alt, "spot_number": spot_number, "consensus_read_count": consensus_read_count_str, "read_quality_A":consensus_read_qual_dict["A"], 
    #           "read_quality_T":consensus_read_qual_dict["T"],"read_quality_C":consensus_read_qual_dict["C"], "read_quality_G":consensus_read_qual_dict["G"]}
    #site_df=pd.DataFrame.from_dict(site_info,orient='index').T
    return site_info


def check_errors(count_dict,ref,threshold=3):
    count_above_threshold = sum(1 for allele in count_dict if allele != ref and count_dict[allele] >= threshold)
    count_nonzero = sum(1 for allele in count_dict if count_dict[allele] > 0)
    pcr_error = 1 if count_above_threshold >= 1 and count_nonzero >= 2 else 0
    pcr_alt=[]
    if pcr_error:
        pcr_alt = [allele for allele in count_dict if allele != ref and count_dict[allele] >= threshold]
        # print("pcr",count_dict)
    
    # lysis_error: require all reads in one UMI must same
    lysis_error = 0
    lysis_alt=[]
    for allele in count_dict:
        if allele != ref and count_dict[allele] > threshold:
            if all(count_dict[other_allele] == 0 for other_allele in count_dict if other_allele != allele):
                lysis_error = 1
                lysis_alt.append(allele)
                # print("lysis",count_dict)
                break

    return pcr_error, lysis_error,pcr_alt,lysis_alt


def combine_UMI_bulk_for_errors(site_barcode_UMI_dict,chrom,pos,ref,threshold):
    #out_df = pd.DataFrame(columns=["chr","pos","ref","alt","spot_number","consensus_read_count","read_quality_A","read_quality_T","read_quality_C","read_quality_G"])
    #spot_number=len(site_barcode_UMI_dict.keys())
    pcr_errors,lysis_errors=0,0
    pcr_alts,lysis_alts=[],[]
    UMI_dp=0
    dp=0
    # lysis_UMI_count=defaultdict(int)
    for barcode in site_barcode_UMI_dict.keys():
        UMI_dp+=len(site_barcode_UMI_dict[barcode].keys())
        for UMI in site_barcode_UMI_dict[barcode]:
            count_dict=site_barcode_UMI_dict[barcode][UMI]["count"]
            dp+=sum(count_dict.values())
            #quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
            pcr_error, lysis_error,pcr_alt,lysis_alt=check_errors(count_dict,ref,threshold)
            pcr_errors += pcr_error; lysis_errors += lysis_error
            pcr_alts += pcr_alt; lysis_alts += lysis_alt

    if lysis_errors >1 or lysis_errors==0:
        lysis_alt_out="."   
    else:
        lysis_alt_out=lysis_alts[0]

    pcr_alt_out=",".join(list(set(pcr_alts)))
    # lysis_alt_out=",".join(list(set(lysis_alts)))
    if pcr_alt_out =="":
        pcr_alt_out="." 
    # if lysis_alt_out =="":
    #     lysis_alt_out="."     

    error_info=[chrom, pos, dp,UMI_dp,ref, pcr_alt_out, pcr_errors, lysis_alt_out,lysis_errors]

    return error_info



### combine UMI for individual data
def UMI_combination_spot(site_barcode_UMI_dict,chrom,pos,ref):
    #out_df = pd.DataFrame(columns=["chr","pos","ref","alt","spot_number","consensus_read_count","read_quality_A","read_quality_T","read_quality_C","read_quality_G"])
    spot_number=len(site_barcode_UMI_dict.keys())
    all_genos=["A","T","C","G"]
    consensus_read_count={}
    consensus_read_quality={}
    for barcode in site_barcode_UMI_dict.keys():
        consensus_read_count[barcode]=defaultdict(int)
        consensus_read_quality[barcode]=defaultdict(dict)
        for UMI in site_barcode_UMI_dict[barcode]:
            count_dict=site_barcode_UMI_dict[barcode][UMI]["count"]
            quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
            phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
            candidate_allele,phred=get_most_candidate_allele(phred_dict,ref)

            consensus_read_count[barcode][candidate_allele]+=1
            if str(phred) not in consensus_read_quality[barcode][candidate_allele].keys():
                consensus_read_quality[barcode][candidate_allele][str(phred)]=0

            consensus_read_quality[barcode][candidate_allele][str(phred)]+=1
    
    # calculate alt
    # alt=",".join([alt_geno for alt_geno in consensus_read_count.keys() if alt_geno != ref])

    # consensus_read_count_list=[]
    # consensus_read_qual_dict=defaultdict(str)
    # for barcode in consensus_read_quality.keys():
    #     for geno in all_genos:
    #         geno_count=0 if geno not in consensus_read_count.keys() else consensus_read_count[geno]
    #         consensus_read_count_list.append(str(geno_count))
            
    #         if geno not in consensus_read_quality.keys():
    #             geno_qual_str="NA" 
    #         else:
    #             geno_qual_list=[str(phred)+":"+str(consensus_read_quality[geno][phred]) for phred in consensus_read_quality[geno].keys()]
    #             geno_qual_str=",".join(geno_qual_list)
    #         consensus_read_qual_dict[geno]=geno_qual_str

    #     consensus_read_count_str=",".join(consensus_read_count_list)

    
    # site_info={"chrom":chrom, "site":pos, "ID":".","ref":ref, "alt":alt, "spot_number": spot_number, "consensus_read_count": consensus_read_count_str, "read_quality_A":consensus_read_qual_dict["A"], 
    #            "read_quality_T":consensus_read_qual_dict["T"],"read_quality_C":consensus_read_qual_dict["C"], "read_quality_G":consensus_read_qual_dict["G"]}
    # site_df=pd.DataFrame.from_dict(site_info,orient='index').T
    return consensus_read_count, consensus_read_quality


def combine_UMI_spot_tidy(bam_dict,barcode_dict,run_type,identifier):
    '''
    identifier:"chr1_1000_A_alt"
    '''
    identifier_list=identifier.split("_")
    #print(identifier_list)
    chr=identifier_list[0]
    pos=int(identifier_list[1])
    ref=identifier_list[2]
    new_list=[]
    all_genos=["A","T","C","G"]
    bam_file=choose_chr_file(bam_dict,chr)
    # print("\n======",chr,bam_file)
    handel_bam=read_bam(bam_file)
    site_barcode_UMI_dict=handle_reads_per_pos_read_count(handel_bam,chr,pos,run_type)
    if barcode_dict=={}:
        filter_site_barcode_UMI_dict=site_barcode_UMI_dict
    else:
        filter_site_barcode_UMI_dict={}
        for key in site_barcode_UMI_dict.keys():
            if key.split(".")[0] in list(barcode_dict.keys()):
                filter_site_barcode_UMI_dict[key]=site_barcode_UMI_dict[key]

    consensus_read_count, consensus_read_quality=UMI_combination_spot(site_barcode_UMI_dict,chr,pos,ref)
    for barcode in consensus_read_count.keys():
        # if ref == "N":
        #     ref=sorted(consensus_read_count[barcode].items(), key = lambda item:item[1], reverse=True)[0][0]

        alt=",".join([alt_geno for alt_geno in consensus_read_count[barcode].keys() if alt_geno != ref])
        if alt=="":
            alt="."
        consensus_read_count_list=[]
        consensus_read_qual_dict=defaultdict(str)
        for geno in all_genos:
            geno_count=0 if geno not in consensus_read_count[barcode].keys() else consensus_read_count[barcode][geno]
            consensus_read_count_list.append(str(geno_count))
            
            if geno not in consensus_read_quality[barcode].keys():
                geno_qual_str="NA" 
            else:
                geno_qual_list=[str(phred)+":"+str(consensus_read_quality[barcode][geno][phred]) for phred in consensus_read_quality[barcode][geno].keys()]
                geno_qual_str=",".join(geno_qual_list)
            consensus_read_qual_dict[geno]=geno_qual_str

        consensus_read_count_str=",".join(consensus_read_count_list)

        barcode_list=[chr,pos,".",ref,alt,barcode,consensus_read_count_str,consensus_read_qual_dict["A"],consensus_read_qual_dict["T"],consensus_read_qual_dict["C"],consensus_read_qual_dict["G"]]
        new_list.append(barcode_list)

    return new_list


def combine_UMI_bulk_tidy(bam_dict,barcode_dict,run_type,identifier):
    '''
    identifier:"chr1_1000_A_alt"
    '''
    identifier_list=identifier.split("_")
    #print(identifier_list)
    chr=identifier_list[0]
    pos=int(identifier_list[1])
    ref=identifier_list[2]
    new_list=[]
    all_genos=["A","T","C","G"]
    bam_file=choose_chr_file(bam_dict,chr)
    # print("\n======",chr,bam_file)
    handel_bam=read_bam(bam_file)
    site_barcode_UMI_dict=handle_reads_per_pos_read_count(handel_bam,chr,pos,run_type)
    if barcode_dict=={}:
        filter_site_barcode_UMI_dict=site_barcode_UMI_dict
    else:
        filter_site_barcode_UMI_dict={}
        for key in site_barcode_UMI_dict.keys():
            if key.split(".")[0] in list(barcode_dict.keys()):
                filter_site_barcode_UMI_dict[key]=site_barcode_UMI_dict[key]

    site_info=combine_UMI_bulk(filter_site_barcode_UMI_dict,chr,pos,ref)

    return site_info


def combine_UMI_bulk_for_errors_tidy(bam_file,barcode_dict,threshold,run_type,identifier):
    '''
    identifier:"chr1_1000_A_alt"
    '''
    identifier_list=identifier.split("_")
    #print(identifier_list)
    chr=identifier_list[0]
    
    pos=int(identifier_list[1])
    ref=identifier_list[2]
    # print("\n======",chr,bam_file)
    handel_bam=read_bam(bam_file)
    site_barcode_UMI_dict,strand=handle_reads_per_pos_read_count_and_strand(handel_bam,chr,pos,run_type)
    if barcode_dict=={}:
        filter_site_barcode_UMI_dict=site_barcode_UMI_dict
    else:
        filter_site_barcode_UMI_dict={}
        for key in site_barcode_UMI_dict.keys():
            if key.split(".")[0] in list(barcode_dict.keys()):
                filter_site_barcode_UMI_dict[key]=site_barcode_UMI_dict[key]

    error_info=combine_UMI_bulk_for_errors(filter_site_barcode_UMI_dict,chr,pos,ref,threshold) + [strand]

    return error_info

# def combine_all(bam_dict,barcode_dict,identifier):
#     spot_info=combine_UMI_spot_tidy(bam_dict,barcode_dict,identifier)
#     ind_info=combine_UMI_bulk_tidy(bam_dict,barcode_dict,identifier)
#     return [spot_info, ind_info]

# @profile   
# import os
from memory_profiler import memory_usage,profile
# @profile(precision=4, stream=open(f"partial_func_pid{os.getpid()}.log", "w"))
def combine_all(bam_dict,barcode_dict,run_type,identifier):
    # mem_now = memory_usage(-1, interval=0.1)[0]
    # print(mem_now)
    # print(identifier,"start!")
    ## maybe we could only calculate spot count, and combine them by cluster or bulk info, based on spot count
    spot_info=combine_UMI_spot_tidy(bam_dict,barcode_dict,run_type,identifier)
    # ind_info=combine_UMI_bulk_tidy(bam_dict,barcode_dict,identifier)
    # print(identifier," done!")
    return [spot_info]


def combine_bulk(bam_dict,barcode_dict,run_type,identifier):
    #spot_info=combine_UMI_spot_tidy(bam_dict,barcode_dict,identifier)
    ind_info=combine_UMI_bulk_tidy(bam_dict,barcode_dict,run_type,identifier)
    return [ind_info]


def combine_bulk_for_errors(bam,barcode_dict,threshold,run_type,mv_chr,identifier):
    #spot_info=combine_UMI_spot_tidy(bam_dict,barcode_dict,identifier)
    # print("type:",run_type)
    
    if identifier.split("_")[0] in mv_chr:
        return []

    else:
        error_info=combine_UMI_bulk_for_errors_tidy(bam,barcode_dict,threshold,run_type,identifier)
    # print(identifier,"done")
        return [error_info]