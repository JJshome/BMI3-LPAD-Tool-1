from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

# 定义标准碱基
STANDARD_BASES = ["A", "C", "G", "T"]

def filter_sequences(input_file, output_file):
    """
    过滤FASTA文件中的非标准序列
    - 去除含有非标准碱基的序列
    - 确保序列长度一致
    - 确保没有插入或缺失

    参数：
    input_file (str): 输入FASTA文件路径
    output_file (str): 输出过滤后FASTA文件路径
    """
    sequences = list(SeqIO.parse(input_file, "fasta"))
    filtered_sequences = []

    print(f"正在过滤 {len(sequences)} 条序列...")
    for record in tqdm(sequences, desc="过滤序列"):
        sequence = str(record.seq).upper()
        
        # 检查是否包含非标准碱基
        if not all(base in STANDARD_BASES for base in sequence):
            continue
        
        # 检查序列长度一致性
        if len(sequence) != len(sequences[0].seq):
            continue

        # 确保没有插入或缺失
        if "-" in sequence:
            continue

        # 添加到过滤后的列表
        filtered_sequences.append(record)
    
    print(f"过滤完成，共保留 {len(filtered_sequences)} 条序列。")

    # 保存过滤后的序列
    SeqIO.write(filtered_sequences, output_file, "fasta")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="过滤FASTA文件中的非标准序列")
    parser.add_argument("--input", required=True, help="输入FASTA文件路径")
    parser.add_argument("--output", required=True, help="输出过滤后FASTA文件路径")
    args = parser.parse_args()

    filter_sequences(args.input, args.output)
