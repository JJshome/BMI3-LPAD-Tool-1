import csv
from calc_features_single import calc_tm, gc_percent, calc_terminal_delta_g, if_secondary_structure
from calc_features_multi import primers_distance, if_dimer
from primer_score import score_primers  # 确保打分函数文件正确导入
from Bio.Seq import Seq


def reverse_complement(sequence):
    """
    获取DNA序列的反向互补。
    :param sequence: 原始DNA序列
    :return: 反向互补序列
    """
    return str(Seq(sequence).reverse_complement())


def parse_primer_info(primer_info):
    """
    从字符串中解析引物信息字典。
    :param primer_info: 字符串格式的引物信息
    :return: 引物信息字典
    """
    return eval(primer_info)  # 确保输入格式正确


def process_input_file(input_file='data/output/Intermediate_file/specific_primer.csv',
                        output_file='data/output/Intermediate_file/score_into_model.csv'):
    """
    处理输入文件，计算每个引物组的打分并保存结果。
    :param input_file: 输入CSV文件路径
    :param output_file: 输出结果文件路径
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        
        # 定义保留字段
        primer_scores = ['gc_score', 'tm_score', 'delta_g_score', 'hairpin_score']
        primers_names = ['F3', 'F2', 'F1', 'B1', 'B2', 'B3']
        fieldnames = (
            ['lamp_id'] +
            [f"{primer}_{score}" for primer in primers_names for score in primer_scores] +
            ['distance_score', 'dimer_score']
        )

        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for idx, row in enumerate(reader):
            # 提取引物信息并调整顺序
            primers = [
                parse_primer_info(row['forward_outer_info'])['sequence'],  # F3
                parse_primer_info(row['forward_middle_info'])['sequence'],  # F2
                parse_primer_info(row['forward_inner_info'])['sequence'],  # F1
                reverse_complement(parse_primer_info(row['reverse_inner_info'])['sequence']),  # B1
                reverse_complement(parse_primer_info(row['reverse_middle_info'])['sequence']),  # B2
                reverse_complement(parse_primer_info(row['reverse_outer_info'])['sequence']),  # B3
            ]

            # 提取起始位置
            start_positions = [
                parse_primer_info(row['forward_outer_info'])['position'],  # F3
                parse_primer_info(row['forward_middle_info'])['position'],  # F2
                parse_primer_info(row['forward_inner_info'])['position'],  # F1
                parse_primer_info(row['reverse_inner_info'])['position'],  # B1
                parse_primer_info(row['reverse_middle_info'])['position'],  # B2
                parse_primer_info(row['reverse_outer_info'])['position'],  # B3
            ]

            # 调用打分函数
            scores = score_primers(primers, start_positions)
            single_primer_scores = scores['single_primer_scores']

            # 构造输出字典，仅保留需要的列
            output_row = {'lamp_id': f'Group_{idx + 1}'}
            for primer, details in zip(primers_names, single_primer_scores):
                for score_key in primer_scores:
                    output_row[f"{primer}_{score_key}"] = details.get(score_key, None)

            # 添加 distance_score 和 dimer_score
            output_row['distance_score'] = scores['distance_score']
            output_row['dimer_score'] = scores['dimer_score']

            # 写入结果
            writer.writerow(output_row)


def save_full_scores(input_file='data/output/Intermediate_file/specific_primer.csv', 
                     full_output_file='data/output/Intermediate_file/full_score_results.csv'):
    """
    保存 score_primers 函数的完整输出到 CSV 文件。
    :param input_file: 输入CSV文件路径
    :param full_output_file: 保存完整分数的输出文件路径
    """
    with open(input_file, 'r') as infile, open(full_output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)

        # 定义完整结果的字段
        fieldnames = ['lamp_id', 'full_scores']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for idx, row in enumerate(reader):
            # 提取引物信息并调整顺序
            primers = [
                parse_primer_info(row['forward_outer_info'])['sequence'],  # F3
                parse_primer_info(row['forward_middle_info'])['sequence'],  # F2
                parse_primer_info(row['forward_inner_info'])['sequence'],  # F1
                reverse_complement(parse_primer_info(row['reverse_inner_info'])['sequence']),  # B1
                reverse_complement(parse_primer_info(row['reverse_middle_info'])['sequence']),  # B2
                reverse_complement(parse_primer_info(row['reverse_outer_info'])['sequence']),  # B3
            ]

            # 提取起始位置
            start_positions = [
                parse_primer_info(row['forward_outer_info'])['position'],  # F3
                parse_primer_info(row['forward_middle_info'])['position'],  # F2
                parse_primer_info(row['forward_inner_info'])['position'],  # F1
                parse_primer_info(row['reverse_inner_info'])['position'],  # B1
                parse_primer_info(row['reverse_middle_info'])['position'],  # B2
                parse_primer_info(row['reverse_outer_info'])['position'],  # B3
            ]

            # 调用打分函数
            scores = score_primers(primers, start_positions)

            # 写入完整结果
            writer.writerow({'lamp_id': f'Group_{idx + 1}', 'full_scores': scores})


if __name__ == "__main__":
    input_csv = 'data/output/Intermediate_file/specific_primer.csv'
    full_output_csv = 'data/output/Intermediate_file/full_score_results.csv'  # 替换为实际完整分数输出路径
    save_full_scores(input_csv, full_output_csv)


if __name__ == "__main__":
    input_csv = 'data/output/Intermediate_file/specific_primer.csv'
    output_csv = 'data/output/Intermediate_file/score_into_model.csv'  # 替换为实际输出文件路径
    process_input_file(input_csv, output_csv)
