# LAMP Primer Auto design(LPAD) Tool

### 1. **设计需求**：
![Pasted image 20241106154754](https://github.com/user-attachments/assets/53370d77-0a43-437c-82fa-0a87d77aef73)

LAMP扩增需要六个引物：前内引物 (FIP)、后内引物 (BIP)、两个外引物 (F3和B3)，以及两个环引物 (LF和LB)。这些引物是用来锁定特定的靶序列区域的。其中，环引物不一定是需要的。本项目主要针对于自动设计前内引物 (FIP)、后内引物 (BIP)、两个外引物 (F3和B3)，并自动化输出优质的F1,F2,F3;B1,B2,B3引物。

### 2. **主要挑战**：

- **复杂的设计限制**：LAMP引物设计比传统的PCR引物复杂，需要考虑多个因素：
    - **特异性**：引物只能结合在目标序列上，不能在其他污染序列上（人类基因组）扩增。
    - **熔解温度 (Tm)**：最佳范围是58-63°C。
    - **GC含量**：理想值是40-60%。
    - **ΔG（自由能）**：需要保持低ΔG，防止引物不稳定
    - **二级结构预测和二聚体**：避免产生发卡结构，二聚体等。
    - **引物间距**：LAMP引物设计要求非常严格的间距和结合模式，比如FIP和BIP要在目标序列上相隔一定的距离，确保高效的环状结构形成，从而实现等温扩增。

### 3. **工具流程**：

- **输入数据**：
    - 工具输入包括靶序列的比对信息（你希望扩增的区域）和背景基因组信息（你不希望扩增的区域）。
- **序列比对和引物生成**：
    - 引物在一定的限制条件下（长度、间距）从靶序列中生成，且严格检查特异性，防止引物与背景（污染）基因组结合。
- **定义限制条件和评分标准**：
    - 对每个候选引物进行约束条件的评估（如Tm、GC含量、自由能），并根据其符合程度赋予分数。
    - 加入二级结构预测和二聚体形成过滤，这将进一步提高评分准确性，剔除潜在问题引物。
- **输出结果**：
    - 工具将输出排好序的引物组合，按分数从高到低列出最佳候选供用户选择。

## 目前思路：

### 算法可行性：

- **TM**:

Tm is estimated using the Nearest-Neighbor method. This method is currently considered to be the
 approximation method that gives the value closest to the actual value. 
 
The calculated Tm is affected by experimental conditions such as the salt concentration and oligo concentration,
 so it is preferred that Tm be calculated under fixed experimental conditions (oligo concentration at 0.1 µM, sodium
 ion concentration at 50 mM, magnesium ion concentration at 4 mM). 
 
The Tm for each region is designed to be about 65°C (64 - 66°C) for F1c and B1c, about 60°C (59 - 61°C) for F2,
 B2, F3, and B3, and about 65°C (64 - 66°C) for the loop primers. 
 
- **GC**%:
  
Primers are designed so that their GC content is between about 40% to 65%.

Primers with GC content between 50% and 60% tend to give relatively good primers. 

- **ΔG**:

The end of the primers serves as the starting point of the DNA synthesis and thus must have certain degree of
stability.   The 3’ ends of F2/B2, F3/B3, and LF/LB and the 5’ end of F1c/B1c are designed so that the free energy
is –4 kcal/ mol or less. The 5’ end of F1c after amplification corresponds to the 3’ end of F1, so that stability is 
important.

- **引物间距**:
  
![image](https://github.com/user-attachments/assets/31048e0f-552f-4ecf-b96f-41a9465b3db9)

（参考Primer Explorer文档）

- **二级结构预测和二聚体**:

考虑序列内部不会产生互补序列(hairpin),以及序列间不会互补(dimer)

hairpin设计：1. 发卡结构的互补区段长度在6-12bp，2. 环状区段的长度范围4-8bp.

dimer设计：检查是否序列间有8-16bp的互补区域.

- **评分系统**:
  
1. 对单个引物的各个特征/以及引物组的特征进行初步评分
2. 使用LAMPPrimerBank数据库数据与nonsense primer数据进行决策树模型训练，确定各特征的权值
3. 引物集(F1-3,B1-3)的总分为各个特征的初评分*权值(第二步获得)

### Benchmark 三步：

1. primer可行性：取任何published文章中使用的primer，投入我们的工具，检查primer的分数是否较好；
2. 工具可行性：取published文章中的扩增序列与选用的primer，将扩增序列投入检查是否结果中有对应primer且分数良好；
3. 各工具对比：取相同扩增序列若干，投入PrimerExplorer5, NEB LAMP and PremierBiosoft，检查各个工具结果是否类似。
