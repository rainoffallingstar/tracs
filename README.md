# tracs

`tracs`（原 `trackplot-rs`）是用 Rust 原生实现 bigWig 取值、整理与**出图**的命令行工具（并兼容替代 `bwtool` 的 `summary/matrix` 子命令）。不依赖 R 或显示服务器。

核心能力：
- `summary`：兼容 `bwtool summary -with-sum -keep-bed -header <bed> <bigwig> <out>`
- `matrix`：兼容 `bwtool matrix -starts/-ends -tiled-averages=<bin> <up:down> <bed> <bigwig> <out>`
- `track-extract`：更高层的一次性提取（多 bigWig + loci/gene + binsize + 可选 gene 模型/ideogram），供 `trackplot.R` 直接消费
- `plot-track`（别名 `plot`）：`track-extract` + **Rust 原生渲染**输出 PDF/SVG（无需 R；默认自动分配配色）
- `profile`：`profile_plot()` 的原生实现（`tracs matrix` 矩阵 → 每个样本一条均值/中位数曲线）
- `heatmap`：`profile_heatmap()` 的原生实现（每个样本一个 panel，按行均值/中位数排序）
- `pca`：`pca_plot()` 的原生实现（`tracs summary` 汇总表 → 样本 PCA 散点图 + 方差解释 scree panel）
- `volcano`：`volcano_plot()` 的原生实现（差异分析结果表 → volcano 图；不重跑 limma，只接受通用的 logFC / p / padj 表）

> **全流程零 R 依赖。** 早期版本通过 `Rscript trackplot.R` 出图；现在布局、绘图、PDF 生成都在 Rust 内完成（`src/plot/`）。构建、测试、CI 都不安装或调用 R，`tracs plot` 也不需要 R、X11 或显示服务器。
> 输出格式由 `--out` 的扩展名决定（`.pdf` 或 `.svg`）。
> `trackplot.R` 仍保留在仓库中，供 R 侧集成（`TRACKTOOLS_TRACKPREP_CMD`）使用，但已不在 `tracs plot` 的调用链上。

## 构建

本环境里 `~/.cargo/config` 把 `crates-io` 替换到了 `rsproxy.cn`，可能无法解析域名。建议用独立的 `CARGO_HOME`：

```bash
cd tracs
CARGO_HOME=/tmp/cargo-home CARGO_REGISTRIES_CRATES_IO_PROTOCOL=sparse cargo build --release
```

产物：`target/release/tracs`

查看帮助：

```bash
./target/release/tracs --help
./target/release/tracs track-extract --help
./target/release/tracs plot --help
```

## 测试（自动化对齐检查）

```bash
cd tracs
CARGO_HOME=/tmp/cargo-home CARGO_REGISTRIES_CRATES_IO_PROTOCOL=sparse cargo test
```

说明：
- CI 会从 GEO `GSE199964` 下载 `GSE199964_RAW.tar` 并解压出若干 `.bigWig/.bw` 作为端到端测试数据源（不提交到 Git）。
- 本地测试默认读取 `localdata/data/GSE199964_RAW/`（可用 `TRACKTOOLS_TESTDATA_DIR` 覆盖），并在固定 loci 上验证 `summary/matrix/track-extract` 的数值正确性。
- 如果系统里额外安装了 `bwtool`（PATH 可找到），测试会再跑一遍 `bwtool summary/matrix` 并逐列对比，作为“旧 bwtool 路径”的对齐校验；没安装则自动跳过该对比用例。
- 如果 `bwtool` 只在 mamba 环境里，可通过设置 `TRACKTOOLS_TEST_BWTOOL_ENV=<env>` 让测试用 `micromamba run -n <env> bwtool ...`（优先）或 `conda run -n <env> bwtool ...` 来完成对齐对比（取决于系统里可用的 runner）。
- 如果你的 micromamba root prefix 不在默认位置，额外设置 `TRACKTOOLS_TEST_MICROMAMBA_ROOT=<root>`（等价于传 `micromamba run -r <root> ...`）。

### 端到端出图测试（`tests/e2e_plot.rs`）

`tests/e2e_plot.rs` 会用真实 bigWig 跑完整的 `tracs plot` 链路（`track-extract` → Rust 原生渲染 → PDF），覆盖：

- 多样本 `--coldata` + `--loci` 出图，并校验 `tracks.tsv` 的 bin 数、样本名与信号非零；
- `--gene`（ENSG）+ **完整** hg19 `ensGene` GTF 出图，校验解析到的染色体/起止/链向与 exon 模型；
- `--gene`（symbol）经 `org.Hs.eg.db` 的 sqlite 本地映射到 ENSG 后出图（只用系统 `sqlite3`，不需要 R）。

依赖与跳过策略（缺任一项即打印说明并跳过，不会让 `cargo test` 失败）：

- `localdata/data/hg19.ensGene.gtf`（完整注释，可用 `TRACKTOOLS_TEST_GTF` 覆盖）；
- symbol 模式还需要本地 orgdb；`TRACKTOOLS_ORGDB_SQLITE` 或 `TRACKTOOLS_GENE_LOOKUP=online` 二者其一。

完整 GTF 可从 UCSC 获取（CI 用同一下载源）：

```bash
curl -L -o /tmp/hg19.ensGene.gtf.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ensGene.gtf.gz
gunzip -f /tmp/hg19.ensGene.gtf.gz && mv /tmp/hg19.ensGene.gtf localdata/data/
```

CI 里会把渲染出的 PDF 与中间 TSV 作为 `e2e-plot-artifacts` 上传，便于失败时排查；本地可设 `TRACKTOOLS_TEST_KEEP_WORKDIR=1` 保留临时工作目录。

准备 `bwtool`（可选，仅用于对齐测试）：

```bash
micromamba create -y -r /tmp/tracktools-mamba -n tracktools-bwtool \
  -c pwwang -c bioconda/label/cf201901 -c conda-forge/label/cf201901 \
  bwtool htslib openssl=1.0.2h
```

运行对齐测试：

```bash
cd tracs
TRACKTOOLS_TEST_BWTOOL_ENV=tracktools-bwtool TRACKTOOLS_TEST_MICROMAMBA_ROOT=/tmp/tracktools-mamba \
  CARGO_HOME=/tmp/cargo-home CARGO_REGISTRIES_CRATES_IO_PROTOCOL=sparse cargo test
```

如果遇到 micromamba cache lock 权限问题，可临时指定：

```bash
export XDG_CACHE_HOME=/tmp/tracktools-xdg-cache
```

## 对接 trackplot.R

有两种方式：

1) 只替换 `bwtool`（`trackplot.R` 仍然在 R 里生成窗口并逐个调用子命令）：

```r
Sys.setenv(TRACKTOOLS_BWTOOL_CMD = "./target/release/tracs")
```

之后 `trackplot.R` 里原本调用 `bwtool summary/matrix` 的地方会自动改用该命令。

2) 让 Rust 接管 `track_extract()` 的“画图前步骤”（窗口生成 + bigWig 取值 + 可选 GTF 查 gene 模型）：

```r
Sys.setenv(TRACKTOOLS_TRACKPREP_CMD = "./target/release/tracs")
```

当 `TRACKTOOLS_TRACKPREP_CMD` 被设置时，`trackplot.R` 的 `track_extract()` 会调用 `track-extract` 并读取 `tracks.tsv/meta.tsv`，不再依赖 conda/bwtool。

兼容性：旧的 `GREATCHIP_BWTOOL_CMD` / `GREATCHIP_TRACKPREP_CMD` 仍然可用，但不再推荐。

说明：
- `track-extract --gene <...>` 默认通过 UCSC `refGene`（Rust 原生 MySQL 客户端）解析基因坐标/外显子（需要网络）。
- 如果提供 `--gtf`，会优先走本地 GTF 做离线查找（推荐在离线环境使用）。
- `track-extract` 默认会从 UCSC 拉取 `cytoBand`（输出 `cytoband.tsv`），用于 `track_plot(show_ideogram=TRUE)` 画 ideogram；不需要 ideogram 时可加 `--no-cytoband` 跳过。

如果你不需要 ideogram，可以在 CLI 里加 `--no-cytoband` 跳过。

`--gene` 输入类型与默认转换（自动校验/归一化）：
- 未指定 `--gtf`：会把输入归一化为 **gene symbol**（支持传 symbol / Entrez / ENSG），再去查 UCSC `refGene.name2`。
- 指定了 `--gtf`：会把输入归一化为 **ENSG**（支持传 symbol / Entrez / ENSG），优先用本地 GTF 做离线查找。
- 本地转换默认使用 `org.Hs.eg.db` 的 sqlite（通过系统 `sqlite3` 读取，**不需要安装 R**）。可用 `TRACKTOOLS_ORGDB_SQLITE=/path/to/org.Hs.eg.sqlite` 覆盖路径。
- 从 Bioconductor 直接取该 sqlite（CI 用的就是这种方式）：
  ```bash
  curl -L -o /tmp/org.Hs.eg.db.tar.gz \
    https://bioconductor.org/packages/3.21/data/annotation/src/contrib/org.Hs.eg.db_3.21.0.tar.gz
  tar -xzf /tmp/org.Hs.eg.db.tar.gz -C /tmp org.Hs.eg.db/inst/extdata/org.Hs.eg.sqlite
  export TRACKTOOLS_ORGDB_SQLITE=/tmp/org.Hs.eg.db/inst/extdata/org.Hs.eg.sqlite
  ```
- 如果没有本地 orgdb，但你希望在线转换，可设置 `TRACKTOOLS_GENE_LOOKUP=online`（需要可联网；使用 mygene.info REST，Rust 内置实现，纯 Rust：reqwest + rustls）。

## Rust 主程序直接出图（推荐）

```bash
./target/release/tracs plot \
  --out out.pdf \
  --gene SLC19A1 --build hg19 --binsize 200 \
  --bigwig localdata/data/GSE199964_RAW/H3K27ac_1.bigWig --sample H3K27ac_1
```

更贴近 `track_plot()` 的一键用法（推荐 `--coldata`）：

```bash
cat > coldata.tsv <<'EOF'
bw_files\tbw_sample_names
localdata/data/GSE199964_RAW/H3K27ac_1.bigWig\tH3K27ac_1
localdata/data/GSE199964_RAW/H3K4me3_1.bigWig\tH3K4me3_1
EOF

./target/release/tracs plot \
  --out out.pdf \
  --loci chr1:158145820-158156686 \
  --binsize 200 \
  --coldata coldata.tsv \
  --col auto \
  --show-ideogram false \
  --draw-gene-track false \
  --group-auto-scale true \
  --show-axis true
```

说明：
- 默认不需要显式指定颜色：`--col auto` 会自动为 tracks 分配离散色盘；也可以传 `--col "#d34,#2980b9,..."` 手动指定。
- 常用 `track_plot()` 参数已映射为 CLI 参数（例如 `--y-max/--y-min`、`--bw-ord`、`--layout-ord`、`--bw-track-height`、`--gene-track-height`、`--cytoband-track-height`、`--regions-bed`、`--boxcol/--boxcolalpha`）。

## profile / heatmap / pca / volcano

这四个子命令分别对应 `trackplot.R` 里的 `profile_plot()`、`profile_heatmap()`、`pca_plot()`、`volcano_plot()`。前三个的输入是 `tracs matrix` 矩阵或 `tracs summary` 汇总表，`volcano` 的输入是差异分析结果表：

```bash
# profile：每个样本一条曲线（输入是 `tracs matrix` 的矩阵，可重复 --matrix）
./target/release/tracs profile \
  --matrix H3K27ac_1.matrix --matrix H3K4me3_1.matrix \
  --sample H3K27ac_1,H3K4me3_1 \
  --up 2500 --down 2500 --stat mean \
  --out profile.pdf

# heatmap：每个样本一个 panel，按行均值/中位数排序
./target/release/tracs heatmap \
  --matrix H3K27ac_1.matrix --matrix H3K4me3_1.matrix \
  --sample H3K27ac_1,H3K4me3_1 \
  --sort-by mean --col-pal Blues \
  --out heatmap.pdf

# pca：样本 PCA 散点图 + 方差解释 scree panel
# 输入是 `tracs summary -with-sum` 的每样本汇总表（可重复 --summary），不是 `matrix` 矩阵
./target/release/tracs pca \
  --summary H3K27ac_1.summary --summary H3K27ac_2.summary --summary H3K4me3_1.summary \
  --sample H3K27ac_1,H3K27ac_2,H3K4me3_1 \
  --color-by H3K27ac,H3K27ac,H3K4me3 \
  --top 1000 --log2 \
  --out pca.pdf
```

说明：
- `pca --summary` 的每行是一个 region，每列是一个样本（即 `extract_summary()` 的 `data` 形状）；`extract_summary()` 只保留每个 `bwtool summary` 输出的 `sum` 列，本实现同样如此。
- `--top` 会先按行标准差降序取前 N 个 region 再做 PCA（与 `pca_plot()` 一致）；`--log2` 对应 `pca_plot(log2 = TRUE)`，先做 `log2(x + --log2-offset)`。
- PCA 的方差解释与得分与 R 的 `prcomp()` 对齐，并有用例锁定（`testdata/pca_r_oracle.tsv`，由 R 4.6.0 生成）。**分量符号不保证一致**：R 文档明确说明 `prcomp()` 的符号是任意的、甚至不同 R 构建之间都可能不同，所以这里使用固定约定（最大载荷取正）；如需匹配某张参考图，可用 `--flip-x/--flip-y`。
- `--work-dir` 会额外写出中间结果便于核对：`heatmap_limits.tsv`、`pca_components.tsv`、`pca_scores.tsv`、`pca_regions.tsv`。

### volcano

```bash
# limma::topTable() 的输出（TSV，含 # contrast: ... 注释头）
./target/release/tracs volcano \
  --results limma_results.tsv \
  --fdr 0.1 \
  --out volcano.pdf

# DESeq2::results() 的输出（CSV，列名不同）
./target/release/tracs volcano \
  --results deseq_results.csv \
  --fdr 0.05 \
  --out volcano_deseq.pdf
```

输入只要是带表头的差异分析结果表，且包含 logFC、p-value、adjusted p-value 三列即可，列名会自动识别：

| 列 | 自动识别的列名（不区分大小写） |
| --- | --- |
| logFC | `logFC`、`log2FoldChange`、`log2FC` |
| p-value | `P.Value`、`pvalue`、`p_value`、`pval` |
| adjusted p-value | `adj.P.Val`、`padj`、`adj.pval`、`p.adjust`、`fdr` |

列名不匹配时可用 `--logfc-col/--p-col/--padj-col` 显式指定。分隔符是 Tab 还是逗号会自动判断。

说明：
- **本实现不重跑 limma**，只画图。`diffpeak()` 依赖 limma 的 empirical Bayes（`lmFit`/`eBayes`）来算出 `P.Value` 与 `adj.P.Val`；这部分没有等价的 Rust 实现，所以在 R 侧（或任何其他工具）算好后把表交给 `tracs volcano` 即可。
- 显著性判定与 `volcano_plot()` 完全一致：`adj.P.Val < fdr`（严格小于），再按 `logFC` 正负分成 up/down。因此 `logFC` 为 0 或缺失的显著 peak 不计入 up/down，图例两个计数之和可能小于显著 peak 总数。
- `--title` 不传时会读取表头注释里的 `contrast`（`topTable()` 保存时会写成 `# contrast: ...`）。
- 与 R 的一处**有意差异**：`volcano_plot()` 用 `xlims = range(res$logFC)`，只要有任何一行的 `logFC` 是 `NA`，R 的 `range()` 就返回 `NA`，随后 `plot()` 直接以 `need finite 'xlim' values` 失败，整张图什么都画不出来。本实现会跳过这些行、正常画出其余 peak，并在 stderr 提示跳过了多少行（`volcano_summary.tsv` 里也有 `skipped` 列）。
- `P.Value` 恰好为 0 时（`-log10(0)` 为无穷）会报错退出，因为 R 在同样输入下也会因 `ylim` 为 `Inf` 而失败；此时需要先过滤或给这些行设一个下限。
- `--work-dir` 会写出 `volcano_summary.tsv`（总数/可绘制数/跳过数/up/down 计数/坐标范围）。

## Release

GitHub Actions 会按日历日期（UTC）发布 Release tag（`YYYYMMDD`），同一天内多次构建会复用同一个 tag，并把不同平台产物作为 assets 上传到同一个 Release。
