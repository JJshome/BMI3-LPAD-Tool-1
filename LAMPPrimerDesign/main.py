import argparse
import logging
import os
from scripts.filter_sequences import filter_sequences
from scripts.calculate_features import calculate_features
from scripts.design_primers import generate_primers
from scripts.score_primers import score_primers
from scripts.align_and_verify import verify_primers
from scripts.visualize_results import visualize_primers
from scripts.generate_report import generate_report

def setup_logging():
    """设置日志配置"""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    return logging.getLogger()

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="LAMP Primer Design Program")
    parser.add_argument("--input", required=True, help="输入目标序列文件路径（FASTA）")
    parser.add_argument("--reference", required=True, help="输入参考基因组文件路径（FASTA）")
    parser.add_argument("--output", required=True, help="输出目录路径")
    parser.add_argument("--config", default="config/params.json", help="配置文件路径")
    parser.add_argument("--visualize", action="store_true", help="生成可视化报告")
    parser.add_argument("--threads", type=int, default=4, help="并行线程数")
    return parser.parse_args()

def main():
    """主程序入口"""
    logger = setup_logging()
    args = parse_args()

    # Step 1: Filter sequences
    filtered_output = os.path.join(args.output, "filtered_sequences.fasta")
    logger.info("Step 1: Filtering sequences...")
    filter_sequences(args.input, filtered_output)
    logger.info(f"Filtered sequences saved to {filtered_output}")

    # Step 2: Calculate sequence features
    features_output = os.path.join(args.output, "sequence_features.csv")
    logger.info("Step 2: Calculating sequence features...")
    calculate_features(filtered_output, features_output)
    logger.info(f"Sequence features saved to {features_output}")

    # Step 3: Generate candidate primers
    candidate_primers_output = os.path.join(args.output, "candidate_primers.csv")
    logger.info("Step 3: Generating candidate primers...")
    generate_primers(features_output, candidate_primers_output)
    logger.info(f"Candidate primers saved to {candidate_primers_output}")

    # Step 4: Score primers
    scored_primers_output = os.path.join(args.output, "scored_primers.csv")
    logger.info("Step 4: Scoring primers...")
    score_primers(candidate_primers_output, features_output, scored_primers_output)
    logger.info(f"Scored primers saved to {scored_primers_output}")

    # Step 5: Verify primers
    verified_primers_output = os.path.join(args.output, "verified_primers.csv")
    verification_log = os.path.join(args.output, "primer_verification.log")
    logger.info("Step 5: Verifying primers...")
    verify_primers(scored_primers_output, args.reference, verified_primers_output, verification_log)
    logger.info(f"Verified primers saved to {verified_primers_output}")
    logger.info(f"Verification log saved to {verification_log}")

    # Step 6: Visualize primers
    visualization_dir = os.path.join(args.output, "visualizations")
    logger.info("Step 6: Visualizing primer distribution and characteristics...")
    visualize_primers(verified_primers_output, visualization_dir)
    logger.info(f"Visualization charts saved to {visualization_dir}")

    # Step 7: Generate report
    report_file = os.path.join(args.output, "LAMP_Primer_Design_Report.pdf")
    logger.info("Step 7: Generating PDF report...")
    generate_report(verified_primers_output, visualization_dir, report_file)
    logger.info(f"Report generated at {report_file}")

    logger.info("Pipeline execution completed.")

if __name__ == "__main__":
    main()
