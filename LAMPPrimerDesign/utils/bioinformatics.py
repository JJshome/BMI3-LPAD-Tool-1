from Bio.SeqUtils import gc_fraction, MeltingTemp as mt

def calculate_gc_content(sequence):
    """计算序列的GC含量"""
    return gc_fraction(sequence)

def calculate_tm(sequence):
    """计算引物的熔解温度"""
    return mt.Tm_NN(sequence)
