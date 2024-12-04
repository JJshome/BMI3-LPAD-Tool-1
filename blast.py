import subprocess
import os

def create_blast_db(fasta_file, db_name):
    """
    使用makeblastdb构建BLAST数据库
    :param fasta_file: 输入的FASTA文件路径
    :param db_name: 输出数据库的名称（不带扩展名）
    """
    try:
        # 构建BLAST数据库
        subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl', '-out', db_name], check=True)
        print(f"BLAST database {db_name} successfully create")
    except subprocess.CalledProcessError as e:
        print(f"Falied to build BLAST database: {e}")

def run_blast(query_sequences, db_name, output_file):
    """
    使用blastn进行比对
    :param query_sequence: 输入的查询序列（单条）
    :param db_name: BLAST数据库名称
    :param output_file: 输出结果的文件
    """
    try:
        # 创建一个临时文件来保存查询序列
        with open("query.fasta", "w") as query_file:
            for seq in query_sequences:
                query_file.write(f">query_sequence\n{seq}")

        # 运行BLAST比对
        subprocess.run(['blastn', '-query', 'query.fasta', '-db', db_name, '-out', output_file, '-outfmt', '6',"-task", "blastn-short", "-word_size", "7", "-evalue", "1"], check=True)
        print(f"Alignment completed, result saved in {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"BLAST alignment failed: {e}")
    finally:
        # 删除临时查询文件
        if os.path.exists("query.fasta"):
            os.remove("query.fasta")

def parse_blast_output(output_file):
    """
    解析BLAST比对结果
    :param output_file: BLAST比对结果文件
    :return: 解析后的比对结果（列表）
    """
    count=0
    try:
        with open(output_file, "r") as f:
            for line in f:
                count+=1
            if count>0:
                return False
            else:
                return True
        print(f"Alignment result successfully parsed.")
    except Exception as e:
        print(f"Failed to parse: {e}")
def main():
    fasta_file = "hg38.fa"  # 输入FASTA文件路径
    query_sequences = ["CCTCCGCCTTCAGCCACTT","TTGCATTAACTGGTTGTCAAACCACCAACAATGGTTCTGGCGGCGG"]  # 输入查询序列
    db_name = "blast_db"  # BLAST数据库名称
    output_file = "blast_output.txt"  # BLAST输出文件路径

    # 1. 构建BLAST数据库
    create_blast_db(fasta_file, db_name)

    # 2. 执行BLAST比对
    run_blast(query_sequences, db_name, output_file)

    # 3. 解析BLAST比对结果
    results = parse_blast_output(output_file)
    
    print(results)

if __name__ == "__main__":
    main()



