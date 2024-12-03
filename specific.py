"""
两种方法：
1. 直接比对target和reference序列，类似primerexplorer那样。最后输出带有比对信息的primers
2. 针对引物评估特异性，调用primer blast等可以做到，每个引物都要去比对，耗时太久
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO

import xml.etree.ElementTree as ET


fasta_file = "oneseq.fasta"  # 文件路径
record = SeqIO.read(fasta_file, "fasta")
result_handle= NCBIWWW.qblast("blastn", "nt", record.seq)
with open("my_blast.xml","w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()




# 解析 BLAST XML 文件
def extract_conserved_sequences(blast_xml_file):
    # 加载 XML 文件
    tree = ET.parse(blast_xml_file)
    root = tree.getroot()
    
    conserved_sequences = []

    
    # 遍历每个 Hit
    for hit in root.findall(".//Hit"):
        
        # 遍历每个 Hsp (High-scoring Segment Pair)
        for hsp in hit.findall(".//Hsp"):
            hsp_conserved=[]
            qseq = hsp.find("Hsp_qseq").text  # 查询序列
            hseq = hsp.find("Hsp_hseq").text  # 目标序列
            midline = hsp.find("Hsp_midline").text  # 对齐信息
            identity = hsp.find("Hsp_identity").text

            if int(identity) <len(record.seq)-60:
                # 提取保守序列（即中线部分为 '|' 的区域）
                conserved_region = []
                for i in range(len(midline)):
                    if midline[i] == "|":  # 保守位置
                        conserved_region.append(qseq[i]) # 添加查询序列的碱基
                    else:
                        conserved_region.append("-")
                
                conserved_sequence = ''.join(conserved_region)
                if conserved_sequence:
                    partial_seq=conserved_sequence.split("-")
                    for part in partial_seq:
                        if len(part)>10:
                            conserved_sequences.append(part)
    
    print(conserved_sequences)
    filtered_conserved_seq=[]
    for s1 in conserved_sequences:
        bool=True
        for s2 in conserved_sequences:
            if s1 in s2 and s1!=s2:
                bool=False
        if bool:
            filtered_conserved_seq.append(s1)
    print(filtered_conserved_seq)
    return filtered_conserved_seq

# 使用自定义的 XML 文件路径
blast_xml_file = "my_blast.xml"
conserved_sequences = extract_conserved_sequences(blast_xml_file)

with open("conserved_sequences.txt","w") as f:
    for seq in conserved_sequences:
        f.write(f">{seq}\n")
        f.write(f"{len(seq)}")
        f.write("\n")

