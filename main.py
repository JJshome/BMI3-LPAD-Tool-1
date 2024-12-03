from primer_score import score_primers
import plot_report

# 假设 primer_list 和 start_pos_list 是你的输入数据
primer_list = ["CTGGTTGTCAAACAACTGG", "TAATAATCTTGGTGGCGTTG", "TACCATAACCAGCTGCTGCG", 
               "TCAAGTGCAAACTCTAGCCGT", "CAGCAGCACCAAGAACTG", "TTCTCTTTCTGGTCAAGGTA"]
start_pos_list = [50, 69, 142, 171, 200, 226]

# 计算打分
scores = score_primers(primer_list, start_pos_list)

# 生成PDF报告
create_pdf_with_features_and_scores(primer_list, start_pos_list, scores)
