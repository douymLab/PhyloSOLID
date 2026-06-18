#!/usr/bin/env python3
import re
import sys

def process_file(input_file, output_file):
    """
    处理输入的annotation文件，生成包含五列的输出文件
    列名: Symbol, Mutation ID, Gene ID, Genome region, Function prediction
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 写入表头
        outfile.write("Symbol\tMutation ID\tGene ID\tGenome region\tFunction prediction\n")
        
        for line in infile:
            line = line.strip()
            if not line:  # 跳过空行
                continue
            
            # 分割每一行，预期有三列
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            # 解析第一列: chr2_3576510_G_A
            mutation_parts = parts[0].split('_')
            if len(mutation_parts) == 4:
                chr_part, pos_part, ref_part, alt_part = mutation_parts
                # 转换为 chr1:39034563,T>A 格式
                mutation_id = f"{chr_part}:{pos_part},{ref_part}>{alt_part}"
            else:
                # 如果不是预期的4部分，保留原样
                mutation_id = parts[0]
            
            # 解析第二列: exonic(synonymous_SNV) 或 ncRNA_exonic
            second_col = parts[1]
            genome_region = ""
            function_pred = "unknown"
            
            # 检查是否有括号
            match = re.search(r'^(.*?)\((.*?)\)$', second_col)
            if match:
                # 有括号的情况
                genome_region = match.group(1)  # 括号前的内容
                function_pred = match.group(2)  # 括号内的内容
            else:
                # 没有括号的情况
                genome_region = second_col
                # function_pred保持为"unknown"
            
            # 第三列: Symbol
            symbol = parts[0]
            
            # 写入输出文件
            outfile.write(f"{symbol}\t{mutation_id}\t{symbol}\t{genome_region}\t{function_pred}\n")

def main():
    # 检查命令行参数
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        print("Example: python script.py Org4S15D63.benchmark_identifier.anno_info.txt output.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        process_file(input_file, output_file)
        print(f"Processing complete! Output written to {output_file}")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
