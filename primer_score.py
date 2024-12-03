import math
from calc_features_single import calc_tm, gc_percent, calc_terminal_delta_g, if_secondary_structure
from calc_features_multi import primers_distance, if_dimer


def score_primers(primer_list, start_pos_list):
    """
    计算6个引物组合的加权评分，并返回可视化报告所需的数据。
    :param primer_list: 6个引物序列的列表
    :param start_pos_list: 6个引物起始位置的列表
    :return: 包含各个引物特征的字典，方便后续绘制图表
    """
    # 权重设置（可以根据具体需求进行调整）
    WEIGHT_GC = 0.2  # GC含量权重
    WEIGHT_TM = 0.3  # Tm值权重
    WEIGHT_DELTA_G = 0.2  # 自由能（Delta G）权重
    WEIGHT_HAIRPIN = 0.3  # 发卡结构权重

    # 初始化返回数据字典
    primer_scores = []
    total_score = 0

    # 定义评分常数
    GC_CONTENT_MIN = 40  # GC含量最低值
    GC_CONTENT_MAX = 65  # GC含量最高值
    GC_CONTENT_OPT = (50, 60)  # GC含量的最优范围
    Tm_F1c_B1c_MIN = 64  # F1c和B1c的Tm最低值
    Tm_F1c_B1c_MAX = 66  # F1c和B1c的Tm最高值
    Tm_F2_B2_F3_B3_LOOP_MIN = 59  # F2, B2, F3, B3, Loop引物的Tm最低值
    Tm_F2_B2_F3_B3_LOOP_MAX = 61  # F2, B2, F3, B3, Loop引物的Tm最高值
    DELTA_G_MAX = -4  # 3'端和5'端的自由能要求（单位：kcal/mol）
    DISTANCE_F2_B2_MIN = 120  # F2到B2的最小距离
    DISTANCE_F2_B2_MAX = 160  # F2到B2的最大距离
    DISTANCE_F2_F1_MIN = 40  # F2到F1的最小距离
    DISTANCE_F2_F1_MAX = 60  # F2到F1的最大距离
    DISTANCE_F2_F3_MAX = 60  # F2到F3的最大距离
    DISTANCE_B2_B1_MIN = 40  # B2到B1的最小距离
    DISTANCE_B2_B1_MAX = 60  # B2到B1的最大距离
    DISTANCE_B2_B3_MAX = 60  # B2到B3的最大距离

    # 计算每个引物的特征
    for i, primer in enumerate(primer_list):
        tm = calc_tm(primer)
        gc = gc_percent(primer)
        terminal_delta_g_5p, terminal_delta_g_3p = calc_terminal_delta_g(primer)
        secondary_structure = if_secondary_structure(primer)

        # 计算引物的单项特征得分
        gc_score = max(0, min(1, 1 - (gc - (GC_CONTENT_MAX + GC_CONTENT_MIN) / 2) / (GC_CONTENT_MAX - GC_CONTENT_MIN)))
        if GC_CONTENT_OPT[0] <= gc <= GC_CONTENT_OPT[1]:
            gc_score += 1
        gc_score = round(gc_score, 2)

        def calculate_tm_score(tm, tm_min, tm_max):
            if tm_min <= tm <= tm_max:
                return 1  # 在范围内打最高分
            elif tm < tm_min - 5:
                return 0  # 超过下限5度，打0分
            elif tm > tm_max + 5:
                return 0  # 超过上限5度，打0分
            else:
                # 计算偏离程度
                if tm < tm_min:
                    # 在下限偏离，构造非线性函数
                    score = 1 - ((tm_min - tm) ** 2 / (5 ** 2))
                else:
                    # 在上限偏离
                    score = 1 - ((tm - tm_max) ** 2 / (5 ** 2))

                return round(max(0, min(1, score)), 2)  # 确保分数在0到1之间

        # 进行tm打分
        if i == 2 or i == 3:  # F1c, B1c
            tm_score = calculate_tm_score(tm, Tm_F1c_B1c_MIN, Tm_F1c_B1c_MAX)
        else:  # F2, B2, F3, B3, Loop
            tm_score = calculate_tm_score(tm, Tm_F2_B2_F3_B3_LOOP_MIN, Tm_F2_B2_F3_B3_LOOP_MAX)
        delta_g_score = 1 if terminal_delta_g_5p <= DELTA_G_MAX and terminal_delta_g_3p <= DELTA_G_MAX else 0
        hairpin_score = 0 if secondary_structure else 1

        # 使用权重来计算加权得分
        weighted_score = round((gc_score * WEIGHT_GC + tm_score * WEIGHT_TM +
                                delta_g_score * WEIGHT_DELTA_G + hairpin_score * WEIGHT_HAIRPIN), 2)

        # 将每个引物的得分保存
        primer_scores.append({
            'sequence': primer,
            'gc': gc,
            'tm': tm,
            'delta_g_5p': terminal_delta_g_5p,
            'delta_g_3p': terminal_delta_g_3p,
            'hairpin': hairpin_score,
            'gc_score': gc_score,
            'tm_score': tm_score,
            'delta_g_score': delta_g_score,
            'hairpin_score': hairpin_score,
            'weighted_score': weighted_score
        })

        # 计算总得分
        total_score += weighted_score

    # 计算引物间的多重特征（距离）
    distances = primers_distance(primer_list, start_pos_list)
    distance_F2_B2_score = round(max(0, min(1, (DISTANCE_F2_B2_MAX - distances['distance_F2_to_B2_end']) / (
                DISTANCE_F2_B2_MAX - DISTANCE_F2_B2_MIN))), 2)
    distance_F2_F1_score = round(max(0, min(1, (DISTANCE_F2_F1_MAX - distances['distance_F2_to_F1']) / (
                DISTANCE_F2_F1_MAX - DISTANCE_F2_F1_MIN))), 2)
    distance_F2_F3_score = round(
        max(0, min(1, (DISTANCE_F2_F3_MAX - distances['distance_F2_to_F3']) / DISTANCE_F2_F3_MAX)), 2)
    distance_B2_B1_score = round(max(0, min(1, (DISTANCE_B2_B1_MAX - distances['distance_B2_to_B1']) / (
                DISTANCE_B2_B1_MAX - DISTANCE_B2_B1_MIN))), 2)
    distance_B2_B3_score = round(
        max(0, min(1, (DISTANCE_B2_B3_MAX - distances['distance_B2_to_B3']) / DISTANCE_B2_B3_MAX)), 2)

    # 计算二聚体评分
    dimer_count = if_dimer(primer_list)
    dimer_score = max(0, 4 - dimer_count)  # 没有二聚体时评分为1，存在二聚体时扣分

    # 计算二聚体状态
    dimer_count = if_dimer(primer_list)

    total_score += distance_F2_B2_score + distance_F2_F1_score + distance_F2_F3_score + distance_B2_B1_score + distance_B2_B3_score + dimer_score

    # 返回所有可视化报告所需的数据
    return {
        'primer_scores': primer_scores,
        'distances': distances,
        'dimer_count': dimer_count,
        'total_score': total_score
    }


if __name__ == "__main__":
    # 测试样例：6个引物和起始位置
    primer_list = ["CTGGTTGTCAAACAACTGG",
                   "TAATAATCTTGGTGGCGTTG",
                   "TACCATAACCAGCTGCTGCG",
                   "TCAAGTGCAAACTCTAGCCGT",
                   "CAGCAGCACCAAGAACTG",
                   "TTCTCTTTCTGGTCAAGGTA"]
    start_pos_list = [50, 69, 142, 171, 200, 226]
    score = score_primers(primer_list, start_pos_list)
    print(f"Total Primer Score: {score}")
