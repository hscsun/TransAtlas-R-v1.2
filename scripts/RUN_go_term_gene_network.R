# 2. 选择要展示的 GO 条目，推荐使用 GO ID
go_terms <- c(
  "GO:0006281",  # DNA repair
  "GO:0000278",  # mitotic cell cycle
  "GO:0006260"   # DNA replication
)

# 3. 指定已有的差异表达结果和本次输出目录
de_file <- paste0(
  "results/NIH3T3/",
  "00_differential_expression/DESeq2_FBS_vs_Control.csv"
)

out_dir <- "results/NIH3T3/go_term_network_example"

stopifnot(
  file.exists("scripts/go_term_gene_network.R"),
  file.exists(de_file)
)

# 4. 设置参数
run_args <- c(
  "scripts/go_term_gene_network.R",
  "--terms", paste(go_terms, collapse = ","),
  "--de", de_file,
  "--species", "mouse",
  "--out", out_dir,
  "--max-genes", "100",
  "--label-top", "60",
  "--up-cutoff", "2.5",
  "--layout", "kk",
  "--width", "10",
  "--height", "9",
  "--title", "",
  "--seed", "1"
)

# 5. 从当前 R 会话启动 Rscript 执行
status <- system2(
  command = file.path(R.home("bin"), "Rscript"),
  args = shQuote(run_args)
)

if (status != 0L) {
  stop("绘图未完成，请查看控制台中的错误信息。")
}

# 6. 查看生成的文件
list.files(out_dir, full.names = TRUE)