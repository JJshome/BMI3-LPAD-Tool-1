import math  

# Dictionary for delta H values  
nn_h = {  
    'AA': 91, 'AC': 65, 'AG': 78, 'AT': 86, 'AN': 80,  
    'CA': 58, 'CC': 110, 'CG': 119, 'CT': 78, 'CN': 91,  
    'GA': 56, 'GC': 111, 'GG': 110, 'GT': 65, 'GN': 85,  
    'TA': 60, 'TC': 56, 'TG': 58, 'TT': 91, 'TN': 66,  
    'NA': 66, 'NC': 85, 'NG': 91, 'NT': 80, 'NN': 80,  
    'aa': 91, 'ac': 65, 'ag': 78, 'at': 86, 'an': 80,  
    'ca': 58, 'cc': 110, 'cg': 119, 'ct': 78, 'cn': 91,  
    'ga': 56, 'gc': 111, 'gg': 110, 'gt': 65, 'gn': 85,  
    'ta': 60, 'tc': 56, 'tg': 58, 'tt': 91, 'tn': 66,  
    'na': 66, 'nc': 85, 'ng': 91, 'nt': 80, 'nn': 80,  
}  

# Dictionary for delta S values  
nn_s = {  
    'AA': 240, 'AC': 173, 'AG': 208, 'AT': 239, 'AN': 215,  
    'CA': 129, 'CC': 266, 'CG': 278, 'CT': 208, 'CN': 220,  
    'GA': 135, 'GC': 267, 'GG': 266, 'GT': 173, 'GN': 210,  
    'TA': 169, 'TC': 135, 'TG': 129, 'TT': 240, 'TN': 168,  
    'NA': 168, 'NC': 210, 'NG': 220, 'NT': 215, 'NN': 203,  
    'aa': 240, 'ac': 173, 'ag': 208, 'at': 239, 'an': 215,  
    'ca': 129, 'cc': 266, 'cg': 278, 'ct': 208, 'cn': 220,  
    'ga': 135, 'gc': 267, 'gg': 266, 'gt': 173, 'gn': 210,  
    'ta': 169, 'tc': 135, 'tg': 129, 'tt': 240, 'tn': 168,  
    'na': 168, 'nc': 210, 'ng': 220, 'nt': 215, 'nn': 203,  
}  

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def gc_percent(primer):  
    """Calculate GC percentage of the primer."""  
    g_count = primer.count('G')  
    c_count = primer.count('C')  
    total_count = len(primer)  
    return round(((g_count + c_count) / total_count) * 100,2)  

def calc_tm_long(dna):  
    """Calculate Tm for long primers (length > 36)."""  
    gc_pct = gc_percent(dna)  
    return round(81.5 + (16.6 * (math.log(50 / 1000.0) / math.log(10))) + (41.0 * (gc_pct / 100)) - (600.0 / len(dna)),2)

def calc_tm_short(primer):  
    """Calculate Tm for short primers (length <= 36)."""  
    dH = 0  
    dS = 108  
    for i in range(len(primer) - 1):  
        pair = primer[i:i+2]  
        dH += nn_h.get(pair, 0)  
        dS += nn_s.get(pair, 0)  
    dH *= -100.0  
    dS *= -0.1  
    tm = (dH / (dS + 1.987 * math.log(100 / 4000000000.0))) - 273.15 + 16.6 * (math.log(50 / 1000.0) / math.log(10))  
    return round(tm)  


def calc_tm(primer):  
    """Calculate Tm based on primer length."""  
    if len(primer) > 36:  
        return calc_tm_long(primer)  
    elif 10 <= len(primer) <= 36:  
        return calc_tm_short(primer)  
    else:  
        return None  
    

def calc_terminal_delta_g(primer, end_length=6):
    """
    计算引物末端区域的自由能(ΔG)，分别评估5'端和3'端的稳定性。
    :param primer: 引物序列
    :param end_length: 末端区域的长度(默认是6个碱基)
    :return: (5'端自由能值, 3'端自由能值)
    """
    # 确保引物长度大于等于 end_length
    if len(primer) < end_length:
        raise ValueError("Primer length must be greater than or equal to the end_length.")
    
    # 计算末端区域序列
    five_prime_seq = primer[:end_length]  # 获取5'端的 end_length 个碱基
    three_prime_seq = primer[-end_length:]  # 获取3'端的 end_length 个碱基

    def calculate_delta_g(sequence):
        """
        计算给定序列的自由能（ΔG）
        :param sequence: 输入序列（5'端或3'端）
        :return: 自由能值
        """
        dH = 0  # 焓变化
        dS = 0  # 熵变化
        
        # 计算序列的 dH 和 dS
        for i in range(len(sequence) - 1):
            pair = sequence[i:i+2]  # 获取相邻的两碱基对
            dH += nn_h.get(pair, 0)  # 查找该对的焓值
            dS += nn_s.get(pair, 0)  # 查找该对的熵值
        
        dH *= -100.0  # 转换单位，单位从 cal/mol 转为 kcal/mol
        dS *= -0.1  # 转换单位，单位从 cal/mol/K 转为 kcal/mol/K
        
        # 计算自由能（ΔG）
        T = 310  # 温度，假设为37°C，即310K
        delta_g = (dH / (dS + 1.987 * math.log(100 / 4000000000.0))) - 273.15 + 16.6 * (math.log(50 / 1000.0) / math.log(10))
        
        return delta_g

    # 分别计算5'端和3'端的ΔG
    five_prime_delta_g = round(calculate_delta_g(five_prime_seq),2)
    three_prime_delta_g = round(calculate_delta_g(three_prime_seq),2)

    return five_prime_delta_g, three_prime_delta_g


def if_secondary_structure(primer: str) -> bool:
    # 定义发卡结构的最小和最大互补区段长度
    min_stem_length = 6
    max_stem_length = 12
    
    # 定义环状区段的长度范围
    min_loop_length = 4
    max_loop_length = 8
    
    # 遍历引物序列，查找自互补区域
    for length in range(min_stem_length, max_stem_length + 1):
        for i in range(len(primer) - length):
            subseq = primer[i:i + length]
            rev_complement_subseq = reverse_complement(subseq)
            
            # 检查剩余的序列中是否包含反向互补序列
            remaining_seq = primer[i + length:]
            index = remaining_seq.find(rev_complement_subseq)
            
            if index != -1:
                # 找到互补序列后，计算第二段序列的起始位置和第一段序列的结束位置的差距
                distance = (i + length + index) - (i + length)
                
                # 检查环状区段的长度是否在限制范围内
                if min_loop_length <= distance <= max_loop_length:
                    return True
    return False
