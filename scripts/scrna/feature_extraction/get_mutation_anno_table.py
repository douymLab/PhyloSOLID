import argparse
import os
import re
import subprocess
import sys

import pandas as pd



def get_intergration_from_identifier_and_one_file(identifier,query_file,index_in_query="0 1 2 3"):
    """
    This function is used to get the intergration result from identifier and one special query file.
    Note: The top 4 columns of query_file must be "chrom pos ref alt"
    
    """
    current_directory = os.path.dirname(os.path.abspath(__file__))
    compare_pl_path=os.path.join(current_directory,"compare_line_for_annovar.pl")
    perl_script=compare_pl_path
    identifier_line="\t".join(identifier.strip().split("_"))
    command=f"echo -e \"{identifier_line}\" |perl {perl_script} - {query_file} {index_in_query}|sort -u"
    # print(command)
    try:
        result=subprocess.check_output(command,text=True,shell=True)
        return result
    except:
        print(f"Something wrong when run the command: {command}")
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
    


def main():
    variant_input=args.prefix+".variant_function"
    exon_variant_input=args.prefix+".exonic_variant_function"
    if not os.path.exists(variant_input) and not os.path.exists(exon_variant_input):
        print(f'Please check {variant_input} and {exon_variant_input}')
        sys.exit()

    mutation_file=open(args.mut,"r")
    mut_list=[i.strip() for i in mutation_file.readlines()]
    out_file=open(args.outfile,"w")
    # print(variant_input,exon_variant_input)

    variant_input_df=pd.read_csv(variant_input,header=None, names=["pos_anno","gene","chrom","pos","pos_2","ref","alt"],sep="\t")
    variant_input_df=variant_input_df.drop_duplicates()
    variant_input_df["pos"]=variant_input_df["pos"].astype(str)
    variant_input_df["identifier"]=variant_input_df["chrom"]+"_"+variant_input_df["pos"]+"_"+variant_input_df["ref"]+"_"+variant_input_df["alt"]
    variant_input_df=variant_input_df.set_index("identifier")

    exon_variant_df=pd.read_csv(exon_variant_input,header=None, names=["linenum","func_anno","trans_info","chrom","pos","pos_2","ref","alt"],sep="\t")
    exon_variant_df=exon_variant_df.drop_duplicates()
    exon_variant_df["pos"]=exon_variant_df["pos"].astype(str)
    exon_variant_df["identifier"]=exon_variant_df["chrom"]+"_"+exon_variant_df["pos"]+"_"+exon_variant_df["ref"]+"_"+exon_variant_df["alt"]
    exon_variant_df=exon_variant_df.set_index("identifier")

    for identifier in mut_list:
        try:
            variant_info=variant_input_df.loc[identifier]
            # print(variant_info)
            anno=variant_info["pos_anno"]
            gene=variant_info["gene"]
            # print(identifier,anno)
        except:
            anno="NONE"
            gene="NONE"

        if anno != "NONE" and anno != "intergenic":
            gene=gene.split("(")[0]
            if anno=="exonic":
                try:
                    exon_info=exon_variant_df.loc[identifier]
                    exon_anno=exon_info["func_anno"]
                    exon_anno=re.sub(" ", "_", exon_anno)
                    anno=anno+"("+exon_anno+")"
                except:
                    exon_anno="NONE"

        out_file.write(f'{identifier}\t{anno}\t{gene}\n')

    out_file.close()

## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--prefix", "-i",required=True,help="the prefix to read input file, such as : a.anno.variant_function; a.anno.exonic_variant_function. The input should be: a.anno ")
parser.add_argument("--mut", "-m",required=True,help="The mutation identifier list")
parser.add_argument("--outfile", "-o",required=True,type=str,help="The path storing the mutinfo")
args = parser.parse_args()
    

if __name__ == '__main__':
    main()

    