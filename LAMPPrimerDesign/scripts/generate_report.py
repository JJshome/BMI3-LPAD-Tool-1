from fpdf import FPDF
import os
import pandas as pd  # 确保导入 pandas 库
from matplotlib import pyplot as plt

class PDFReport(FPDF):
    def header(self):
        self.set_font("Arial", "B", 12)
        self.cell(0, 10, "LAMP Primer Design Report", 0, 1, "C")

    def footer(self):
        self.set_y(-15)
        self.set_font("Arial", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}", 0, 0, "C")

def generate_report(primer_file, output_dir, output_file):
    """
    生成PDF报告，包含引物特性和可视化图表
    参数：
        primer_file (str): 验证后的引物文件路径
        output_dir (str): 可视化图表存储目录
        output_file (str): 输出PDF报告路径
    """
    # 加载引物数据
    primers = pd.read_csv(primer_file)  # 使用pandas读取数据

    pdf = PDFReport()
    pdf.add_page()

    # 添加概要
    pdf.set_font("Arial", size=12)
    pdf.cell(0, 10, "Summary of Primer Design", ln=1)
    pdf.ln(5)
    pdf.set_font("Arial", size=10)
    pdf.cell(0, 10, f"Total Primers Designed: {len(primers)}", ln=1)
    pdf.cell(0, 10, f"Average GC Content: {primers['GC_Content'].mean():.2f}", ln=1)
    pdf.cell(0, 10, f"Average Tm: {primers['Tm'].mean():.2f}", ln=1)
    pdf.ln(10)

    # 插入可视化图表
    for chart in ["gc_content_distribution.png", "tm_distribution.png", "primer_distribution.png"]:
        chart_path = os.path.join(output_dir, chart)
        if os.path.exists(chart_path):
            pdf.cell(0, 10, chart.replace("_", " ").capitalize(), ln=1)
            pdf.image(chart_path, x=10, y=pdf.get_y(), w=180)
            pdf.ln(65)

    # 保存PDF
    pdf.output(output_file)
    print(f"PDF报告已生成：{output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="生成LAMP引物设计报告")
    parser.add_argument("--input", required=True, help="输入验证后的引物文件路径")
    parser.add_argument("--chart_dir", required=True, help="可视化图表存储目录")
    parser.add_argument("--output", required=True, help="输出PDF报告路径")
    args = parser.parse_args()

    generate_report(args.input, args.chart_dir, args.output)
