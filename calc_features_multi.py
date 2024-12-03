def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def primers_distance(primer_list, start_pos_list):
    # 引物长度
    lengths = [len(primer) for primer in primer_list]
    
    # 计算所需的距离
    distance_F2_to_B2_end = (start_pos_list[4] + lengths[4]) - (start_pos_list[1])  # F2到B2的末端
    distance_F2_to_F1 = start_pos_list[2] - start_pos_list[1]  # F2的5'端到F1的5'端
    distance_F2_to_F3 = start_pos_list[1] - (start_pos_list[0] + lengths[0])  # F2的5'端到F3的末端
    distance_B2_to_B1 = (start_pos_list[4] + lengths[4]) - (start_pos_list[3] + lengths[3])  # B2的末端到B1的末端
    distance_B2_to_B3 = start_pos_list[5] - (start_pos_list[4]+lengths[4])  # B2的末端到B3的5'端
    
    return {
        "distance_F2_to_B2_end": distance_F2_to_B2_end,
        "distance_F2_to_F1": distance_F2_to_F1,
        "distance_F2_to_F3": distance_F2_to_F3,
        "distance_B2_to_B1": distance_B2_to_B1,
        "distance_B2_to_B3": distance_B2_to_B3
    }


def if_dimer(primer_list):
    # 定义互补区段长度的最小和最大阈值
    min_dimer_length = 8
    max_dimer_length = 16
    
    # 用来统计符合条件的dimer对的数量
    dimer_count = 0
    
    # 遍历每对引物（包括引物自己与自己）
    for i in range(len(primer_list)):
        for j in range(i, len(primer_list)):
            primer1 = primer_list[i]
            primer2 = primer_list[j]
            
            # 获取primer2的反向互补序列（用于检查互补性）
            rev_complement_primer2 = reverse_complement(primer2)
            
            # 检查是否存在互补区段
            for length in range(min_dimer_length, max_dimer_length + 1):
                # 查找primer1和rev_complement_primer2之间的互补区段
                for k in range(len(primer1) - length + 1):
                    subseq1 = primer1[k:k + length]
                    if subseq1 in rev_complement_primer2:
                        dimer_count += 1  # 找到符合条件的互补区段，增加计数
                        break  # 一旦找到一个互补区段就停止查找该引物对
    
    return dimer_count  # 返回统计的二聚体数量