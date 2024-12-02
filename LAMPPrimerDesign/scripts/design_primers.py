import pandas as pd
from tqdm import tqdm

def generate_primers(feature_file, output_file, primer_length=20):
    """
    根据序列特征生成候选引物
    参数：
        feature_file (str): 输入的特征矩阵文件路径（CSV格式）
        output_file (str): 输出候选引物文件路径（CSV格式）
        primer_length (int): 引物长度（默认20）
    """
    # 读取特征矩阵
    features = pd.read_csv(feature_file)
    primers = []

    print(f"正在生成引物，每个序列生成多个候选引物...")
    for _, row in tqdm(features.iterrows(), total=features.shape[0], desc="生成引物"):
        sequence = row["ID"]
        seq_length = row["Length"]
        full_sequence = row["ID"]  # 假设输入包含完整序列

        # 根据序列长度生成引物
        for start in range(seq_length - primer_length + 1):
            primer_seq = full_sequence[start:start + primer_length]
            primers.append({
                "Sequence_ID": sequence,
                "Primer_Sequence": primer_seq,
                "Start_Position": start,
                "Length": primer_length
            })

    # 转换为DataFrame
    df_primers = pd.DataFrame(primers)
    df_primers.to_csv(output_file, index=False)
    print(f"候选引物已保存到 {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="生成候选引物")
    parser.add_argument("--input", required=True, help="输入特征矩阵文件路径（CSV格式）")
    parser.add_argument("--output", required=True, help="输出候选引物文件路径（CSV格式）")
    parser.add_argument("--length", type=int, default=20, help="引物长度（默认20）")
    args = parser.parse_args()

    generate_primers(args.input, args.output, args.length)
