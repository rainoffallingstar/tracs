# localdata layout

本仓库把**本地数据**与**分析输出**放在 `localdata/`，并在 Git 中忽略该目录（见 `.gitignore`）。

- `localdata/data/`：输入数据与参考注释（例如 `GSE199964_RAW/`、`hg19.ensGene.gtf`、`gtf.bed`、`test.RDS`）
- `localdata/results/`：分析/作图输出（例如 `tracktools_out/`、各类 PDF、deeptools/macs3 输出等）
- `localdata/scripts/`：项目外的临时脚本/一次性脚本
- `localdata/rstudio/`：RStudio 的本地状态（`.Rproj.user`、`_Rproj.user`、`.Rhistory`）
- `localdata/bin/`：下载的第三方可执行文件（如有）

