# Faecalibacterium 项目最终整理版

本目录是源项目 `/scratch/p312334/project/8--Faecalibacterium` 的非破坏性整理与审计结果。源文件没有被移动或覆盖；小型结果复制到这里，大型原始矩阵和数据库通过符号链接保留出处。

## 最终主数据链

```text
FP_Summary_V4.csv：48,395 genomes
        ↓ CheckM2 completeness > 80 且 contamination < 5
17,104 project-defined high-quality genomes
        ↓ skani triangle + average linkage，95% ANI cut
84 skani_95_average groups = 84 species
        ↓ 每组 1 个代表基因组
84 representatives
```

最终物种键固定为 `skani_95_average`。`Haoran_species` 只有 83 个唯一字符串，是因为 group 19 和 group 23 被赋予了同一个 GTDB 标准化展示名；两组 ANI 为 94.20%，仍应作为两个独立 species。另一个旧列 `skani_group` 也有 84 组，但有 3 个 genome 与 `skani_95_average` 的归属不同，不能混用。

17,104 在本项目中可称为“高质量基因组”，但方法部分必须同时写清项目阈值；它不等同于完整的 MIMAG high-quality MAG 定义。

## V4 中各 resource 的基因组数

| Resource | V4 genomes | 高质量 genomes |
|---|---:|---:|
| 00SPIRE | 31,476 | 9,543 |
| UHGG | 7,205 | 4,278 |
| DAG3 | 6,485 | 1,267 |
| Lifeline_Next | 1,909 | 1,002 |
| GTDB | 904 | 703 |
| NCBI | 167 | 167 |
| 1000IBD | 141 | 66 |
| Africa | 98 | 69 |
| UMCG | 10 | 9 |
| **合计** | **48,395** | **17,104** |

## 13 项分析的整理状态

| 目录 | 状态 | 可用结论 |
|---|---|---|
| [01_genome_collection](01_genome_collection/README.md) | `FINAL` | V4 resource 数量及高质量数量已重新计数、作图 |
| [02_species_rarefaction_qc](02_species_rarefaction_qc/README.md) | `FINAL` | 84 species、rarefaction、CheckM2、84 representatives 已验证 |
| [03_ani_umap_phylogeny](03_ani_umap_phylogeny/README.md) | `FINAL_WITH_DECLARED_EXCLUSIONS` | 17,104-genome ANI dendrogram；UMAP 嵌入 17,076，明确排除 28 个低连通 genome |
| [04_animal_phylogeny](04_animal_phylogeny/README.md) | `REDONE_REPRODUCIBLE` | 17,104-tip GTDB-Tk tree 已按上传的 971-tip 表精确剪枝；人和动物 genome 位于共享物种谱系，支持跨宿主共享/生态联系；单棵树不能证明传播方向 |
| [05_gene_sharing](05_gene_sharing/README.md) | `FINAL_25_GROUP_SUBSET` | 3,968-genome pangenome 的 25 个标签已一对一重标到 final-84 的 25 个 group；可展示 gene sharing，但不是 84 组全覆盖且不能称 HGT |
| [05a_functional_annotation](05a_functional_annotation/README.md) | `REDONE_NUCLEOTIDE_CAMPAIGN` | 源数据有两套 gapseq 运行——蛋白输入走 blastp 打 Prokka 蛋白组，核酸输入走 **tblastn 打 contig**——共用同一个绝对 bitscore 200 阈值而打分不可互换，**不可混用**。推荐核酸口径（4,498 genomes / 19 groups）。物种功能分化与精氨酸–腐胺路线对调（lineage marker，非宿主效应）两口径一致；**宿主效应两口径 strict 结果零重叠、不可报告**，2026-08-19 的 β-氧化结论已撤回；表型 split 只报 prevalence 差值不给 p 值 |
| [06_representatives_sylph_benchmark](06_representatives_sylph_benchmark/README.md) | `LOD_VALIDATED_RANDOM1000_RUNNING` | 85-run 固定梯度/LoD benchmark 已验证；按原方法重建的 final-84、随机 5-group、1,000-replicate 主 benchmark 正在运行（主数组 30649953；task 34 重试 30652262；最终评估 30652810） |
| [06a_early_life_mother_infant](06a_early_life_mother_infant/README.md) | `EXPANDED_84_CONTEXT` | 4,577 个 LLNEXT profile、53 个实际物化 species rows；总体检出率由 2 周 3.02% 上升到 12 月 75.31%，总体及物种轨迹均带 Wilson 95% CI，并显示 6 个全局 prevalence≤5% 的重复低频早期检出物种 |
| [07_coabundance_network](07_coabundance_network/README.md) | `FINAL_DESCRIPTIVE` | 报告“在 nominally significant cohorts 中”的正/负方向比例，不解释为直接互作 |
| [08_geography](08_geography/README.md) | `FINAL_DESCRIPTIVE` | 37,666 个去重 stool samples；continent/country 均为统一 muted-red 气泡矩阵，旧非 slide 图已删除 |
| [09_phenotype_meta_analysis](09_phenotype_meta_analysis/README.md) | `REDONE_FULL_CHAIN_VALIDATED` | 修复 CAD 多层级分母后重跑；abundance 与 prevalence 分别做 endpoint-wide BH-FDR，另以固定规则归纳 23 个跨端点同向 broad axes |
| [09a_bile_acid_pangenome](09a_bile_acid_pangenome/README.md) | `NEW_VALIDATED_25_GROUP_SUBSET` | 胆汁酸候选基因呈物种/菌株偏好；*F. duncaniae* 富集 3β-HSDH 模块，但未检出完整 canonical bai 或 UDCA 7-epimerization 通路；基因分布与 meta 效应的 154 项相关均未通过 FDR |
| [10_blood_group](10_blood_group/README.md) | `ORIGINAL_FIGURE_EXACTLY_REPRODUCED` | 正式 A/B 图与旧脚本在数据变换、分组、检验、颜色和版式上一致，tests TSV byte-identical；协变量调整两部分模型保留为补充 |
| [11_intraspecies_association](11_intraspecies_association/README.md) | `REDONE_NULL_IL6` | 709 个双预处理 ID 已显式统一选 kneaddata；IL6 未复现，旧 succinate/炎症机制解释不成立 |

## 建议先看

1. [完整技术审计报告](00_audit_index/AUDIT_REPORT.md)
2. [机器可读状态表](00_audit_index/STATUS_MATRIX.tsv)
3. [尚需处理事项](00_audit_index/ACTION_REGISTER.tsv)
4. [跨产物验证](00_audit_index/tables/REDO_VALIDATION.tsv)
5. [11/11 statistical fallacy scan](00_audit_index/FALLACY_SCAN.md)
6. [最终图件 QA 与推荐入口](00_audit_index/FIGURE_QA.md)
7. [文件清单](00_audit_index/MANIFEST.tsv)
8. [原始 slide 与配色说明](00_audit_index/reference/SLIDE_STYLE_GUIDE.md)

## 结果表述边界

- Animal tree 支持动物肠道是 *Faecalibacterium* 基因组多样性的重要栖息地，并显示人–动物共享物种级谱系；不能仅凭该树断言近期传播、方向或 reservoir。
- Gene sharing 是共享基因库的描述性 proxy，不等于 HGT，也不能单独证明物种分化的因果机制。
- Co-abundance 的 “always positive/negative” 仅指达到 nominal `p<0.05` 的 cohort 中方向一致。
- 正式 fixed-bin benchmark 中，Sylph 在全部 84 组的 `0.01%`、`0.1%`、`5%` bins 均检出，84-group 零填充 precision 为 1、false positives 为 0；`0.0001%` 与 `0.001%` 均未检出。`0.001%` 约对应中位 195 条 150-bp reads，`0.01%` 约对应中位 1,950 条，因此当前模拟支持的经验阈值是约 `0.01%`，不是 `0.001%`。MetaPhlAn4 名称 crosswalk 仅保守唯一解析 21/84 组，序列经验映射仅得到 17/84 个唯一 SGB，且 18 组发生 SGB collision；但在可唯一解析组中，本设计所测的 `0.01%`、`0.1%`、`5%` 档检出和定量总体良好。该结果说明定制 Sylph workflow 对扩展物种集合覆盖更完整，不等于证明 MetaPhlAn4 软件算法失效。
- 当前 IL6 重分析为阴性；不能继续使用旧 IL6–succinate–inflammation 机制链。
- **05a 的宿主功能差异不可报告**：gapseq 蛋白口径与核酸口径的 strict 结果零重叠且方向不同，无法与注释流程选择分离；2026-08-19 基于混用两种口径报告的「脂肪酸 β-氧化在动物基因组中缺失」已撤回（注释假象）。
- 05a 的两套 gapseq 运行来自数据方（`p318243`），不是本项目生成；核酸口径用 tblastn 直接打 contig，不依赖 Prokka 基因预测，蛋白口径完全建立在 Prokka call 之上。
- 05a 的功能–表型关联：跨 13 个物种的 163,625 次相关检验为阴性且零假设校准良好（小 p 值数量不超过随机期望）。「同一 phenotype 上方向相反的物种功能是否不同」改为只报告 prevalence 差值、不给 p 值——每侧最多 4 个物种，置换 p 的下限（0.167–0.25）已高于任何 FDR 阈值。差值本身很大（最大 96.9 个百分点），但 BMI/DCA/Age 的正向侧都只有 1 个物种，构造上等同于该物种的 repertoire 特征。
- 05a 中纯动物 group 与纯人 group 的「独有功能」绝大多数是 lineage marker：它们在混合宿主物种内部并不区分宿主。引用该表必须带 `verdict` 列，并写明是哪个 gapseq campaign。

目录策略见 [ORGANIZATION_POLICY.md](00_audit_index/ORGANIZATION_POLICY.md)。

当前推荐输出一律位于各模块 `artifacts/final/`；被替代结果只放 `artifacts/archive/`。这一规则也适用于后续更新。
