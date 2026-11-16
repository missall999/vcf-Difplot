# 示例输出：sampleA vs sampleB 差异位点分析

本目录是一套完整的实测输出，用来展示 `vcf_difplot.R` 在真实数据上的效果。
样品名已做脱敏处理（`sampleA` / `sampleB` / `sampleC`）。

## 数据来源

| 项 | 说明 |
|---|---|
| 参考样本 | `sampleA` 的 gVCF（848 MB） |
| 联合分型结果 | `sampleA-sampleB-sampleC.{snp,indel}_clean.vcf.gz` |
| 参考基因组 | 水稻 323 v7.0（`Osativa_323_v7.0.fa`） |
| 测试文件 | `ALL.{snp,indel}.vcf.gz`（sampleB + sampleC） |

**说明**：`sampleA` 的 gVCF 全部记录为 `<NON_REF>`（无实际变异），不能直接与测试文件合并。
因此使用 `sampleA + sampleB + sampleC` 三样本联合分型（GenotypeGVCFs）的成品，
避免重复跑数小时的联合分型。

联合分型 VCF 的样本顺序：`sampleA, sampleB, sampleC`。

## 处理流程

```bash
# 1. GATK 导出制表符格式（3 样本）
gatk VariantsToTable \
  -R Osativa_323_v7.0.fa \
  -V sampleA-sampleB-sampleC.indel_clean.vcf.gz \
  -F CHROM -F POS -GF GT \
  -O indel.table

# 2. 绘图：sampleA 为基准，sampleB 为比较，双端纯合检测
Rscript vcf_difplot.R \
  -i indel.table \
  -b sampleA -c sampleB \
  -o indel_sampleA_vs_sampleB_homo.pdf \
  --baseHetcheck --copHetcheck \
  --preset nature \
  --density \
  --output_table indel.tab \
  --table_min_len_diff 1

# SNP 同理（不加 --table_min_len_diff，因为 SNP 的 LEN_DIFF 恒为 0）
```

## 结果

### SNP

| 指标 | 数值 |
|---|---|
| 输入记录 | 3,924,314 |
| 去除缺失/畸形基因型 | 3,527,734 |
| 双端纯合过滤后 | 2,768,453 |
| **差异位点** | **531,983** |
| 非差异位点 | 2,236,470 |
| 涉及染色体 | 14（Chr1–12, ChrUn, ChrSy） |
| 运行耗时 | 7.01 s / 1.25 GB |

### INDEL

| 指标 | 数值 |
|---|---|
| 输入记录 | 655,910 |
| **差异位点** | **112,704** |
| 涉及染色体 | 14 |

### INDEL 差异位点按染色体分布

| 染色体 | 差异 indel 数 |
|---|---|
| Chr1 | 10,960 |
| Chr2 | 11,633 |
| Chr3 | **27,008** |
| Chr4 | 10,303 |
| Chr5 | 11,357 |
| Chr6 | 3,874 |
| Chr7 | 6,897 |
| Chr8 | 9,628 |
| Chr9 | 2,508 |
| Chr10 | **14,819** |
| Chr11 | 1,571 |
| Chr12 | 1,890 |
| ChrSy | 87 |
| ChrUn | 169 |
| **合计** | **112,704** |

Chr3（27,008）和 Chr10（14,819）差异最集中。

## 输出文件

| 文件 | 说明 |
|---|---|
| `snp_diff_positions.tsv` | SNP 差异位点 + 基因型，531,983 行，5 列（SNP 的 LEN_DIFF 恒为 0，故未过滤） |
| `indel_sampleA_vs_sampleB_homo.pdf` | indel 差异图 |
| `indel_sampleA_vs_sampleB_homo_density.pdf` | indel 密度图 |
| `snp_sampleA_vs_sampleB_homo.pdf` | SNP 差异图 |
| `snp_sampleA_vs_sampleB_homo_density.pdf` | SNP 密度图 |
| `*.png` | 上述 PDF 的 PNG 预览 |

### snp_diff_positions.tsv 格式

```tsv
CHROM	POS	sampleA.GT	sampleB.GT	LEN_DIFF
Chr1	247146	A	C	0
Chr1	272773	T	C	0
Chr1	286819	C	T	0
Chr1	286822	A	G	0
Chr1	286825	G	A	0
```

| 列 | 说明 |
|---|---|
| `CHROM` | 染色体 |
| `POS` | 位置（整数，非科学计数法） |
| `sampleA.GT` | 基准样本基因型（纯合时只显示单等位基因） |
| `sampleB.GT` | 比较样本基因型（同上） |
| `LEN_DIFF` | 两样本等位基因长度的绝对差值（bp） |

**基因型格式**：纯合基因型折叠为单等位基因——`G/G` → `G`，`GTTTTT|GTTTTT` → `GTTTTT`。
杂合/多等位基因型（未启用纯合过滤时）保留完整形式如 `A/ATTT`。

**LEN_DIFF**：`abs(nchar(sampleA 等位基因) - nchar(sampleB 等位基因))`。
本文件为 SNP 位点，两条等位基因等长，因此 `LEN_DIFF` 恒为 0。

indel 侧使用 `--table_min_len_diff 1`，**已剔除 853 条等长替换**（LEN_DIFF = 0），
最终 111,851 条，LEN_DIFF 范围 1–272 bp（该明细表未随示例保留）。

**LEN_DIFF 分布**（过滤后 111,851 条）：

| LEN_DIFF | 数量 | | LEN_DIFF | 数量 |
|---|---|---|---|---|
| 1 | 48,686 | | 5 | 3,276 |
| 2 | 18,489 | | 6 | 3,334 |
| 3 | 9,415 | | 7 | 2,075 |
| 4 | 6,130 | | 8 | 2,185 |

绝大多数（48,686 / 44%）是 1bp 差异。最大差异 272 bp。

## 数据校验

- `snp_diff_positions.tsv`：531,983 行，全部唯一，已排序 ✅
- SNP 与 INDEL 位点重叠数：**0**（两类变异互斥，符合预期）✅
- 基因型折叠：纯合基因型已折叠为单等位基因，无 `/` 或 `|` 重复 ✅
- LEN_DIFF 抽样核对：`G` vs `GTTTTT` = 5、`GA` vs `G` = 1、`A` vs `AGAC` = 3、`GAT` vs `G` = 2、58bp vs 20bp = 38，全部正确 ✅
- 纯合性：抽样 20 万条基因型，homo=200000 / het=0，符合双端纯合过滤 ✅

## 复现命令

```bash
# 在 ms 上
cd ~/vcfdiff_run
source ~/miniconda3/etc/profile.d/conda.sh && conda activate Gatk

gatk VariantsToTable -R /home/user/WorkDir/TQassemble/Osativa_323_v7.0.fa \
  -V ~/data/TMP_VCF/sampleA-sampleB-sampleC.indel_clean.vcf.gz \
  -F CHROM -F POS -GF GT -O indel.table

Rscript vcf_difplot.R -i indel.table -b sampleA -c sampleB \
  -o indel_sampleA_vs_sampleB_homo.pdf \
  --baseHetcheck --copHetcheck --preset nature --density \
  --output_table indel.tab --table_min_len_diff 1
# 说明：indel.tab 是 --output_table 生成的位点明细表，本示例目录未保留该文件。
```
