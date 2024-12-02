import pandas as pd

def score_primers(primer_file, feature_file, output_file, tm_threshold=50.0, gc_threshold=0.4):
    """
    根据熔解温度 (Tm) 和 GC 含量对候选引物进行评分
    参数：
        primer_file (str): 输入候选引物文件路径（CSV格式）
        feature_file (str): 输入特征矩阵文件路径（CSV格式）
        output_file (str): 输出评分引物文件路径（CSV格式）
        tm_threshold (float): Tm 值阈值（默认50.0）
        gc_threshold (float): GC 含量阈值（默认0.4）
    """
    primers = pd.read_csv(primer_file)
    features = pd.read_csv(feature_file)

    # 合并特征信息
    merged = pd.merge(primers, features, left_on="Sequence_ID", right_on="ID")

    # 对每个引物评分
    print("正在评分引物...")
    merged["Score"] = (
        (merged["Tm"] >= tm_threshold).astype(int) +
        (merged["GC_Content"] >= gc_threshold).astype(int)
    )

    # 筛选高分引物
    scored_primers = merged[merged["Score"] > 1]
    scored_primers.to_csv(output_file, index=False)
    print(f"评分后的引物已保存到 {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="对候选引物进行评分")
    parser.add_argument("--input", required=True, help="输入候选引物文件路径（CSV格式）")
    parser.add_argument("--features", required=True, help="输入特征矩阵文件路径（CSV格式）")
    parser.add_argument("--output", required=True, help="输出评分引物文件路径（CSV格式）")
    parser.add_argument("--tm_threshold", type=float, default=50.0, help="Tm 值阈值（默认50.0）")
    parser.add_argument("--gc_threshold", type=float, default=0.4, help="GC 含量阈值（默认0.4）")
    args = parser.parse_args()

    score_primers(args.input, args.features, args.output, args.tm_threshold, args.gc_threshold)
