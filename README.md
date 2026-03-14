# tracktools

Rust 实现的最小 `bwtool` 替代，用于本仓库的 `trackplot.R`（只覆盖用到的子命令）：

- `summary`：对应 `bwtool summary -with-sum -keep-bed -header <bed> <bigwig> <out>`
- `matrix`：对应 `bwtool matrix -starts/-ends -tiled-averages=<bin> <up:down> <bed> <bigwig> <out>`
- `track-extract`：更高层的一次性提取（多 bigWig + loci/gene + binsize），供 `trackplot.R` 直接调用
- `plot-track`：Rust 作为主程序，内部调用 `track-extract`，然后调用 `Rscript` 输出 PDF

## 构建

本环境里 `~/.cargo/config` 把 `crates-io` 替换到了 `rsproxy.cn`，可能无法解析域名。建议用独立的 `CARGO_HOME`：

```bash
cd .
CARGO_HOME=/tmp/cargo-home CARGO_REGISTRIES_CRATES_IO_PROTOCOL=sparse cargo build --release
```

产物：`target/release/tracktools`

## 测试（自动化对齐检查）

```bash
cd .
CARGO_HOME=/tmp/cargo-home CARGO_REGISTRIES_CRATES_IO_PROTOCOL=sparse cargo test
```

说明：
- 测试会用仓库自带的 `localdata/data/GSE199964_RAW/H3K27ac_1.bigWig` 在固定 loci 上验证 `summary/matrix/track-extract` 的数值正确性。
- 如果系统里额外安装了 `bwtool`（PATH 可找到），测试会再跑一遍 `bwtool summary/matrix` 并逐列对比，作为“旧 bwtool 路径”的对齐校验；没安装则自动跳过该对比用例。
- 如果 `bwtool` 只在 mamba 环境里，可通过设置 `TRACKTOOLS_TEST_BWTOOL_ENV=<env>` 让测试用 `micromamba run -n <env> bwtool ...`（优先）或 `conda run -n <env> bwtool ...` 来完成对齐对比（取决于系统里可用的 runner）。
- 如果你的 micromamba root prefix 不在默认位置，额外设置 `TRACKTOOLS_TEST_MICROMAMBA_ROOT=<root>`（等价于传 `micromamba run -r <root> ...`）。

准备 `bwtool`（可选，仅用于对齐测试）：

```bash
micromamba create -y -r /tmp/tracktools-mamba -n tracktools-bwtool \
  -c pwwang -c bioconda/label/cf201901 -c conda-forge/label/cf201901 \
  bwtool htslib openssl=1.0.2h
```

运行对齐测试：

```bash
cd .
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
Sys.setenv(TRACKTOOLS_BWTOOL_CMD = "./target/release/tracktools")
```

之后 `trackplot.R` 里原本调用 `bwtool summary/matrix` 的地方会自动改用该命令。

2) 让 Rust 接管 `track_extract()` 的“画图前步骤”（窗口生成 + bigWig 取值 + 可选 GTF 查 gene 模型）：

```r
Sys.setenv(TRACKTOOLS_TRACKPREP_CMD = "./target/release/tracktools")
```

当 `TRACKTOOLS_TRACKPREP_CMD` 被设置时，`trackplot.R` 的 `track_extract()` 会调用 `track-extract` 并读取 `tracks.tsv/meta.tsv`，不再依赖 conda/bwtool。

兼容性：旧的 `GREATCHIP_BWTOOL_CMD` / `GREATCHIP_TRACKPREP_CMD` 仍然可用，但不再推荐。

说明：
- `track-extract --gene <SYMBOL>` 默认通过 UCSC `refGene`（Rust 原生 MySQL 客户端）解析基因坐标/外显子（需要网络）。
- 如果你传的是 Ensembl 基因 ID（例如 `ENSG...`），可以额外提供 `--gtf` 走本地 GTF，不需要网络。
- `track-extract` 默认也会从 UCSC 拉取 `cytoBand`（输出 `cytoband.tsv`），用于 `track_plot(show_ideogram=TRUE)` 画 ideogram。

如果你不需要 ideogram，可以在 CLI 里加 `--no-cytoband` 跳过。

`--gene` 输入类型与默认转换（自动校验/归一化）：
- 未指定 `--gtf`：会把输入归一化为 **gene symbol**（支持传 symbol / Entrez / ENSG），再去查 UCSC `refGene.name2`。
- 指定了 `--gtf`：会把输入归一化为 **ENSG**（支持传 symbol / Entrez / ENSG），优先用本地 GTF 做离线查找。
- 本地转换默认使用 `org.Hs.eg.db` 的 sqlite（通过系统 `sqlite3` 读取）。可用 `TRACKTOOLS_ORGDB_SQLITE=/path/to/org.Hs.eg.sqlite` 覆盖路径。
- 如果没有本地 orgdb，但你希望在线转换，可设置 `TRACKTOOLS_GENE_LOOKUP=online`（需要可联网；使用 mygene.info REST，Rust 内置实现）。

## Rust 主程序直接出图（推荐）

```bash
./target/release/tracktools plot \
  --trackplot-r trackplot.R \
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

./target/release/tracktools plot \
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

常用 `track_plot()` 参数已映射为 CLI 参数（例如 `--y-max/--y-min`、`--bw-ord`、`--layout-ord`、`--bw-track-height`、`--gene-track-height`、`--cytoband-track-height`、`--regions-bed`、`--boxcol/--boxcolalpha`）。
