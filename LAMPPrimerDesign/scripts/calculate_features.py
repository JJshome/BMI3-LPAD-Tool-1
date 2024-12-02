from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
from tqdm import tqdm

def calculate_shannon_entropy(sequence):
    """
    计算序列的Shannon熵，用于评估序列复杂性
    参数：
        sequence (str): 输入序列
    返回：
        entropy (float): Shannon熵值
    """
    from collections import Counter
    import math

    # 统计碱基频率
    counts = Counter(sequence)
    total_bases = len(sequence)
    entropy = 0.0

    for base, count in counts.items():
        probability = count / total_bases
        entropy -= probability * math.log2(probability)

    return entropy

def calculate_gc_content(sequence):
    """
    计算给定DNA序列的GC含量
    参数：
        sequence (str): 输入的DNA序列（如 'ATGC'）
    返回：
        gc_content (float): 序列中的GC含量，范围是0到1之间
    """
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    total_bases = len(sequence)

    # 计算GC含量比例
    gc_content = (g_count + c_count) / total_bases
    return gc_content

def calculate_features(input_file, output_file):
    """
    计算FASTA文件中每个序列的特征（Tm值、GC含量、Shannon熵）
    参数：
        input_file (str): 输入FASTA文件路径
        output_file (str): 输出特征矩阵文件路径（CSV格式）
    """
    from Bio import SeqIO

    # 读取序列
    records = list(SeqIO.parse(input_file, "fasta"))
    features = []

    print(f"正在计算 {len(records)} 条序列的特征...")
    for record in tqdm(records, desc="计算特征"):
        sequence = str(record.seq).upper()

        # 计算特征
        tm_value = mt.Tm_NN(sequence)  # 熔解温度
        gc_content_value = calculate_gc_content(sequence)  # 使用自定义函数计算GC含量
        entropy = calculate_shannon_entropy(sequence)  # Shannon熵

        features.append({
            "ID": record.id,
            "Length": len(sequence),
            "Tm": tm_value,
            "GC_Content": gc_content_value,
            "Shannon_Entropy": entropy
        })

    # 转换为DataFrame
    df = pd.DataFrame(features)
    df.to_csv(output_file, index=False)
    print(f"特征矩阵已保存到 {output_file}")



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="计算序列特征")
    parser.add_argument("--input", required=True, help="输入FASTA文件路径")
    parser.add_argument("--output", required=True, help="输出特征矩阵文件路径（CSV格式）")
    args = parser.parse_args()

    calculate_features(args.input, args.output)
