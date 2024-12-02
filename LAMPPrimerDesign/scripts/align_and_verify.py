from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
from tqdm import tqdm

def verify_primers(primer_file, reference_file, output_file, log_file, max_mismatches=2):
    """
    验证引物的特异性
    参数：
        primer_file (str): 输入评分引物文件路径（CSV格式）
        reference_file (str): 输入参考基因组文件路径（FASTA格式）
        output_file (str): 输出验证后的引物文件路径（CSV格式）
        log_file (str): 验证日志文件路径
        max_mismatches (int): 最大允许的错配数（默认2）
    """
    # 加载引物数据
    primers = pd.read_csv(primer_file)

    # 加载参考基因组序列
    reference = next(SeqIO.parse(reference_file, "fasta")).seq

    # 设置全局比对参数
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.gap_score = -1

    verified_primers = []
    log_entries = []

    print("正在验证引物的特异性...")
    for _, primer in tqdm(primers.iterrows(), total=primers.shape[0], desc="验证引物"):
        primer_seq = primer["Primer_Sequence"]
        
        # 强制将primer_seq转换为字符串类型
        primer_seq = str(primer_seq)  # 确保 primer_seq 是字符串

        # 比对引物到参考基因组
        alignment = aligner.align(reference, primer_seq)[0]
        mismatches = sum(
            1 for a, b in zip(alignment.query, alignment.target) if a != b and b != "-"
        )

        # 检查错配数是否在允许范围内
        if mismatches <= max_mismatches:
            verified_primers.append(primer)
            log_entries.append(f"PASS: {primer_seq} -> {mismatches} mismatches")
        else:
            log_entries.append(f"FAIL: {primer_seq} -> {mismatches} mismatches")

    # 保存验证后的引物
    verified_primers_df = pd.DataFrame(verified_primers)
    verified_primers_df.to_csv(output_file, index=False)

    # 保存日志
    with open(log_file, "w") as log:
        log.write("\n".join(log_entries))

    print(f"验证完成：验证结果已保存到 {output_file}，日志已保存到 {log_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="验证引物的特异性")
    parser.add_argument("--input", required=True, help="输入评分引物文件路径（CSV格式）")
    parser.add_argument("--reference", required=True, help="输入参考基因组文件路径（FASTA格式）")
    parser.add_argument("--output", required=True, help="输出验证后的引物文件路径（CSV格式）")
    parser.add_argument("--log", required=True, help="验证日志文件路径")
    parser.add_argument("--max_mismatches", type=int, default=2, help="最大允许的错配数（默认2）")
    args = parser.parse_args()

    verify_primers(args.input, args.reference, args.output, args.log, args.max_mismatches)
