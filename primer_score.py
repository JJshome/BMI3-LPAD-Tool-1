import math
from calc_features_single import calc_tm, gc_percent, calc_terminal_delta_g, if_secondary_structure
from calc_features_multi import primers_distance, if_dimer

def score_primers(primer_list, start_pos_list):
    """
    计算6个引物组合的评分，评分考虑单引物特征和引物间的多重特征。
    :param primer_list: 6个引物序列的列表
    :param start_pos_list: 6个引物起始位置的列表
    :return: 总分数
    """
    # 初始化分数
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
    DISTANCE_F2_F1_MIN = 40   # F2到F1的最小距离
    DISTANCE_F2_F1_MAX = 60   # F2到F1的最大距离
    DISTANCE_F2_F3_MAX = 60   # F2到F3的最大距离

    # 计算单引物特征（Tm、GC含量、自由能等）
    for i, primer in enumerate(primer_list):
        tm = calc_tm(primer)
        gc = gc_percent(primer)
        terminal_delta_g = calc_terminal_delta_g(primer)
        secondary_structure = if_secondary_structure(primer)

        # 打印各个引物的特征
        print(f"Primer {i+1}: {primer}")
        print(f"  GC Content: {gc}%")
        print(f"  Tm: {tm}°C")
        print(f"  Terminal Delta G: {terminal_delta_g} kcal/mol")
        print(f"  Hairpin: {'Yes' if secondary_structure else 'No'}")
        print()

        # 评分：GC含量
        gc_score = max(0, min(1, (gc - GC_CONTENT_MIN) / (GC_CONTENT_MAX - GC_CONTENT_MIN)))
        if GC_CONTENT_OPT[0] <= gc <= GC_CONTENT_OPT[1]:
            gc_score += 1  # 如果GC含量在最优范围内，加分
        
        # 评分：Tm
        if i == 0 or i == 3:  # F1c, B1c
            if Tm_F1c_B1c_MIN <= tm <= Tm_F1c_B1c_MAX:
                tm_score = 1
            else:
                tm_score = 0
        else:  # F2, B2, F3, B3, Loop
            if Tm_F2_B2_F3_B3_LOOP_MIN <= tm <= Tm_F2_B2_F3_B3_LOOP_MAX:
                tm_score = 1
            else:
                tm_score = 0

        # 评分：3'端自由能
        if terminal_delta_g <= DELTA_G_MAX:
            delta_g_score = 1
        else:
            delta_g_score = 0

        # 评分：发卡结构
        hairpin_score = 0 if secondary_structure else 1

        # 累加单引物评分
        total_score += gc_score + tm_score + delta_g_score + hairpin_score

    # 计算引物间的多重特征（距离和二聚体）
    distances = primers_distance(primer_list, start_pos_list)
    distance_F2_B2_score = max(0, min(1, (DISTANCE_F2_B2_MAX - distances['distance_F2_to_B2_end']) / (DISTANCE_F2_B2_MAX - DISTANCE_F2_B2_MIN)))
    distance_F2_F1_score = max(0, min(1, (DISTANCE_F2_F1_MAX - distances['distance_F2_to_F1']) / (DISTANCE_F2_F1_MAX - DISTANCE_F2_F1_MIN)))
    distance_F2_F3_score = max(0, min(1, (DISTANCE_F2_F3_MAX - distances['distance_F2_to_F3']) / DISTANCE_F2_F3_MAX))

    # 计算二聚体评分
    dimer_count = if_dimer(primer_list)
    dimer_score = max(0, 1 - dimer_count)  # 没有二聚体时评分为1，存在二聚体时扣分

    # 累加多重特征评分
    total_score += distance_F2_B2_score + distance_F2_F1_score + distance_F2_F3_score + dimer_score

    return total_score


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
