import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os  # 导入 os 模块来检查并创建目录

def visualize_primers(primer_file, output_dir):
    """
    可视化引物的分布和特性
    参数：
        primer_file (str): 输入验证后的引物文件路径（CSV格式）
        output_dir (str): 输出图表的目录
    """
    # 确保输出目录存在，如果不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 读取引物数据
    primers = pd.read_csv(primer_file)

    # 绘制GC含量分布
    plt.figure(figsize=(10, 6))
    sns.histplot(primers["GC_Content"], kde=True, bins=20, color="skyblue")
    plt.title("GC Content Distribution")
    plt.xlabel("GC Content")
    plt.ylabel("Frequency")
    plt.savefig(f"{output_dir}/gc_content_distribution.png")
    plt.close()

    # 绘制Tm值分布
    plt.figure(figsize=(10, 6))
    sns.histplot(primers["Tm"], kde=True, bins=20, color="orange")
    plt.title("Tm Distribution")
    plt.xlabel("Melting Temperature (Tm)")
    plt.ylabel("Frequency")
    plt.savefig(f"{output_dir}/tm_distribution.png")
    plt.close()

    # 绘制引物分布图
    plt.figure(figsize=(12, 6))
    sns.scatterplot(x="Start_Position", y="GC_Content", hue="Score", data=primers, palette="viridis", size="Tm")
    plt.title("Primer Distribution")
    plt.xlabel("Start Position")
    plt.ylabel("GC Content")
    plt.legend(title="Score", loc="upper right")
    plt.savefig(f"{output_dir}/primer_distribution.png")
    plt.close()

    print(f"可视化图表已保存到 {output_dir}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="可视化引物分布和特性")
    parser.add_argument("--input", required=True, help="输入验证后的引物文件路径（CSV格式）")
    parser.add_argument("--output_dir", required=True, help="输出图表的目录")
    args = parser.parse_args()

    visualize_primers(args.input, args.output_dir)
