import subprocess
import pandas as pd
import ast

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

def run_blast(db_name, output_file):
    """
    使用blastn进行比对
    :param query_sequence: 输入的查询序列（单条）
    :param db_name: BLAST数据库名称
    :param output_file: 输出结果的文件
    """
    try:
        # 运行BLAST比对
        #subprocess.run(['blastn', '-query', 'query.fasta', '-db', db_name, '-out', output_file, '-outfmt', '6',"-num_threads", "8"], check=True)
        subprocess.run(['blastn', '-query', 'query.fasta', '-db', db_name, '-out', output_file, '-outfmt', '6',"-task", "blastn-short", "-word_size", "7", "-evalue", "1","-num_threads", "8"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"BLAST alignment failed: {e}")


def parse_blast_output(output_file):
    """
    解析BLAST比对结果
    :param output_file: BLAST比对结果文件
    :return: 解析后的比对结果（列表）
    """
    row_record=[]
    col_record=[]
    try:
        with open(output_file, "r") as f:
            for line in f:
                columns = line.split()
                row_now=int(columns[0].split(",")[0])
                primer_len=int(columns[0].split(",")[1])
                if row_now not in row_record and int(float(columns[2]))==100 and int(columns[3])==primer_len:
                    row_record.append(row_now)
        return row_record
    except Exception as e:
        print(f"Failed to parse: {e}")



def convert_to_dict(x):
    return ast.literal_eval(x) 



def main():
    fasta_file = "hg38.fa"  # 输入FASTA文件路径
    db_name = "blast_db"  # BLAST数据库名称
    output_file = "blast_output.txt"  # BLAST输出文件路径

    df=pd.read_csv("result.csv")
    columns = df.columns
    info_columns = [col for col in columns if 'info' in col]
    row,col=df.shape
    df = pd.read_csv('result.csv', converters={col: convert_to_dict for col in info_columns})

    # 1. 构建BLAST数据库
    create_blast_db(fasta_file, db_name)

    # 2. 执行BLAST比对
    with open("query.fasta", "w") as query_file:
        for i in range(row):
            for j in range(1,9):
                dict=df.iat[i,j]
                seq=dict["sequence"]
                length=dict["length"]
                query_file.write(f">{i},{length},{j}\n{seq}\n")
    run_blast(db_name,output_file)
    row_sub = parse_blast_output(output_file)
    print(row_sub)
    result=df.drop(row_sub) 
    result.to_csv("specific_primer.csv",index=False)

if __name__ == "__main__":
    main()