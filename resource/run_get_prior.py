from functools import partial
import multiprocessing
import argparse
import datetime
import os
# from pybedtools import BedTool

from module.read_file import get_freq_for_mutation_simple_from_annovar, get_freq_for_mutation_simple_from_gnomad, handel_and_split_vcf_for_bed, handle_prior_bed
from utils import check_dir, get_chrom_list_from_file, intersect_pos_gnomad, intersect_pos_gnomad_only_2file


def get_prior(args):
    print("start:",datetime.datetime.now())
    # tmp_dir=os.path.join(args.outdir,"tmp")
    # vcf2bed_split_dir=os.path.join(args.outdir,"tmp","split_vcf")
    # short_gnomad_dir=os.path.join(args.outdir,"tmp","short_gnomad")
    # check_dir(tmp_dir)
    vcf2bed_split_dir=os.path.join(args.outdir,"split_vcf")
    short_gnomad_dir=os.path.join(args.outdir,"short_gnomad")
    check_dir(vcf2bed_split_dir);check_dir(short_gnomad_dir)
    
    # in this step, the vcf file have been turned to bed, and split by chrom
    if args.posfile and not args.vcffile:
        pos_bed_file_dict=handel_and_split_vcf_for_bed(args.posfile,vcf2bed_split_dir)

    elif args.vcffile and not args.posfile:
        pos_bed_file_dict=handel_and_split_vcf_for_bed(args.vcffile,vcf2bed_split_dir)
    else:
        print("The posfile and vcffile is exclusive, please check your input!")
        # sys.exit()

    # print(pos_bed_file_dict)
    print("handel_vcf:", datetime.datetime.now())
    # mutation_identifier_list=handle_vcf(args.vcffile)
    #in_name=args.vcffile.split("/")[-1].split(".")[0]
    in_name=args.sample
    #thread=get_process_cores(args.thread)
    thread=args.thread

    print(f"The following steps will be runned in {thread} thread.")
    # new_list=[]
    # for mutation in mutation_identifier_list: 
    #     new_list.append("_".join(mutation))
    out_prior=open(os.path.join(args.outdir, in_name + args.outname), "w")

    if args.gnomad and not args.annovar:
        gnomad_dict=get_chrom_list_from_file(args.gnomad)
        short_gnomad_dict=intersect_pos_gnomad(pos_bed_file_dict, gnomad_dict, short_gnomad_dir) 
        print("intersect_gnomad_pos:",datetime.datetime.now())
        # print(short_gnomad_dict)
        partial_func=partial(get_freq_for_mutation_simple_from_gnomad,short_gnomad_dict)
    elif args.annovar and not args.gnomad:
        gnomad_dict=get_chrom_list_from_file(args.annovar)
        print("intersect_gnomad_pos:",datetime.datetime.now())
        short_gnomad_dict=intersect_pos_gnomad(pos_bed_file_dict, gnomad_dict, short_gnomad_dir)
        # print(short_gnomad_dict)
        partial_func=partial(get_freq_for_mutation_simple_from_annovar,short_gnomad_dict)
    else:
        print("The gnomad file and annovar file have been redendant, please check your input")
        # sys.exit()

    # gnomad_file=os.path.join(args.gnomad,"gnomad.genomes.v3.1.2.sites."+chrom+".vcf.bgz")
    # partial_func=partial(get_freq_for_mutation_simple_from_gnomad,gnomad_dict)
    print("before_multiprocessing:",datetime.datetime.now())

    for chrom in pos_bed_file_dict.keys():
        bed_file=pos_bed_file_dict[chrom]
        new_list=handle_prior_bed(bed_file)
        print("read_vcf:"+chrom,datetime.datetime.now())
        with multiprocessing.Pool(thread) as pool:
            result=pool.map(partial_func, new_list, chunksize=1000) # type: ignore # per result is (chr, pos, ref, fA, fT, fC,fG)
        

        print("finish_miltiprocessing:"+chrom,datetime.datetime.now())
        for line in result:
            out="\t".join([str(i) for i in line])
            out_prior.write(f'{out}\n')
        print("finish_writing_file:"+chrom,datetime.datetime.now())

    out_prior.close()



def get_prior_2file(args):
    print("start:",datetime.datetime.now())
    short_gnomad_dir=os.path.join(args.outdir,"tmp")
    check_dir(args.outdir); check_dir(short_gnomad_dir)

    bed_file=args.bedfile
    in_name=args.sample

    thread=args.thread

    print(f"The following steps will be runned in {thread} thread.")

    out_prior=open(os.path.join(args.outdir, in_name + args.outname), "w")
    index_name=in_name+args.outname

    if args.gnomad and not args.annovar:
        gnomad_file=args.gnomad
        short_gnomad_file=intersect_pos_gnomad_only_2file(bed_file, gnomad_file, index_name,short_gnomad_dir) 
        print("intersect_gnomad_pos:",datetime.datetime.now())
        # print(short_gnomad_dict)
        partial_func=partial(get_freq_for_mutation_simple_from_gnomad,short_gnomad_file)
    elif args.annovar and not args.gnomad:
        gnomad_file=args.annovar
        print("intersect_gnomad_pos:",datetime.datetime.now())
        short_gnomad_dict=intersect_pos_gnomad_only_2file(bed_file, gnomad_file,index_name, short_gnomad_dir) 
        # print(short_gnomad_dict)
        partial_func=partial(get_freq_for_mutation_simple_from_annovar,short_gnomad_dict)
    else:
        print("The gnomad file and annovar file have been redendant, please check your input")
        # sys.exit()

    # gnomad_file=os.path.join(args.gnomad,"gnomad.genomes.v3.1.2.sites."+chrom+".vcf.bgz")
    # partial_func=partial(get_freq_for_mutation_simple_from_gnomad,gnomad_dict)
    print("before_multiprocessing:",datetime.datetime.now())

    new_list=handle_prior_bed(bed_file)
    # with multiprocessing.Pool(thread) as pool:
    #     result=pool.map(partial_func, new_list, chunksize=3) # type: ignore # per result is (chr, pos, ref, fA, fT, fC,fG)
    
    result=[]
    for line in new_list:
        per_result=partial_func(line)
        result.append(per_result)

    # print("finish_miltiprocessing:",datetime.datetime.now())
    for line in result:
        out="\t".join([str(i) for i in line])
        out_prior.write(f'{out}\n')
    print("finish_writing_file:",datetime.datetime.now())

    out_prior.close()


if __name__ == '__main__':
    ## parameters
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')
    
    # get prior
    parser_prior = subparsers.add_parser('prior', help='get prior from gnomad file')
    parser_prior.add_argument("--sample","-s", required=True, help="sample name")
    parser_prior.add_argument("--posfile",required=False, help="input mutation information POS file; chr pos ref alt")
    parser_prior.add_argument("--vcffile",required=False, help="input mutation information VCF file; chr pos ref alt")
    parser_prior.add_argument("--outdir",required=False, default="./", help="output dir")
    parser_prior.add_argument("--outname",required=False, default=".prior.out", help="output suffix")
    parser_prior.add_argument("--gnomad",required=False, help="gnomad list file containing all gnomad file, or one specific gnomad file")
    parser_prior.add_argument("--annovar",required=False, help="gnomad list file from annovar containing all gnomad file, or one specific file")
    parser_prior.add_argument("--thread",required=False, default=2,type=int, help="the thread you want to use, please make sure thread is equal to the cpu numbers you used" )
    parser_prior.set_defaults(func=get_prior)


    # get prior
    parser_prior_s = subparsers.add_parser('splitprior', help='get prior from gnomad file')
    parser_prior_s.add_argument("--sample","-s", required=True, help="sample name")
    parser_prior_s.add_argument("--bedfile",required=False, help="input mutation information VCF file; chr pos-1 pos ref alt")
    parser_prior_s.add_argument("--outdir",required=False, default="./", help="output dir")
    parser_prior_s.add_argument("--outname",required=False, default=".prior.out", help="output index, if split, just use split index, each will be added")
    parser_prior_s.add_argument("--gnomad",required=False, help="gnomad list file containing all gnomad file, or one specific gnomad file")
    parser_prior_s.add_argument("--annovar",required=False, help="gnomad list file from annovar containing all gnomad file, or one specific file")
    parser_prior_s.add_argument("--thread",required=False, default=2,type=int, help="the thread you want to use, please make sure thread is equal to the cpu numbers you used" )
    parser_prior_s.set_defaults(func=get_prior_2file)

    args = parser.parse_args()
    args.func(args)