import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from primer_score import score_primers


# 测试样例：6个引物和起始位置
primer_list = ["CTGGTTGTCAAACAACTGG",
                "TAATAATCTTGGTGGCGTTG", 
                "TACCATAACCAGCTGCTGCG", 
                "TCAAGTGCAAACTCTAGCCGT", 
                "CAGCAGCACCAAGAACTG", 
                "TTCTCTTTCTGGTCAAGGTA"]
start_pos_list = [50, 69, 142, 171, 200, 226]

primer_data = score_primers(primer_list,start_pos_list)['primer_scores']

# 提取数据
primers = ['F3','F2','F1c','B1c','B2','B3']
gcs = [primer['gc'] for primer in primer_data]
tms = [primer['tm'] for primer in primer_data]
delta_g_5p = [primer['delta_g_5p'] for primer in primer_data]
delta_g_3p = [primer['delta_g_3p'] for primer in primer_data]
weighted_scores = [primer['weighted_score'] for primer in primer_data]

# 设置绘图风格
sns.set(style="whitegrid")

# 创建 PDF 文件来保存报告
with PdfPages('primer_report.pdf') as pdf:
    
    # Tm 的 barplot
    plt.figure(figsize=(8, 6))
    sns.barplot(x=primers, y=tms, palette='Blues_d')
    plt.title('Tm of Primers')
    plt.ylabel('Tm (°C)')
    plt.xlabel('Primer')
    pdf.savefig()
    plt.close()
    
    # GC 的 barplot
    plt.figure(figsize=(8, 6))
    sns.barplot(x=primers, y=gcs, palette='Greens_d')
    plt.title('GC Content of Primers')
    plt.ylabel('GC (%)')
    plt.xlabel('Primer')
    pdf.savefig()
    plt.close()
    
    # ΔG 5' 的 barplot
    plt.figure(figsize=(8, 6))
    sns.barplot(x=primers, y=delta_g_5p, palette='RdBu_r')
    plt.title('ΔG at 5\' End of Primers')
    plt.ylabel('ΔG (kcal/mol)')
    plt.xlabel('Primer')
    pdf.savefig()
    plt.close()

    # ΔG 3' 的 barplot
    plt.figure(figsize=(8, 6))
    sns.barplot(x=primers, y=delta_g_3p, palette='RdBu_r')
    plt.title('ΔG at 3\' End of Primers')
    plt.ylabel('ΔG (kcal/mol)')
    plt.xlabel('Primer')
    pdf.savefig()
    plt.close()

    # Weighted Score 的 barplot
    plt.figure(figsize=(8, 6))
    sns.barplot(x=primers, y=weighted_scores, palette='Purples_d')
    plt.title('Weighted Scores of Primers')
    plt.ylabel('Score')
    plt.xlabel('Primer')
    pdf.savefig()
    plt.close()

    # 创建一个表格展示所有引物的特征
    import pandas as pd
    primer_df = pd.DataFrame(primer_data)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=primer_df.values, colLabels=primer_df.columns, loc='center')
    plt.title('Primer Features Table')
    pdf.savefig()
    plt.close()

print("Report has been generated and saved as 'primer_report.pdf'")
