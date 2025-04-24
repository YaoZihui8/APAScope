apalyzer_quant <- function(out_dir) {
  ### APAlyzer Quant ###

  library(Rsubread)
  library(APAlyzer)

  ## 3UTR Quant ##
  refUTRraw <- read.table(paste0(out_dir,"/updUTRraw_hg19.tsv"), header = TRUE)

  # Quant each sample
  bamdir <- "XXX"
  flsall <- dir(bamdir, pattern = ".Aligned.sortedByCoord.out.bam")
  flsall <- paste0(bamdir, flsall)  # 完整路径
  names(flsall) <- dir(bamdir, pattern = ".Aligned.sortedByCoord.out.bam")
  # 设置输出目录
  output_dir <- out_dir

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  for (file in flsall) {
    identifier <- basename(file)
    UTRdbraw <- REF3UTR(refUTRraw)
    DFUTRraw <- PASEXP_3UTR(UTRdbraw, file, Strandtype = "forward")
    output <- file.path(output_dir, paste0(identifier, ".csv"))
    write.csv(DFUTRraw, file = output, row.names = FALSE, col.names = TRUE)
    cat("Processed:", identifier, "\n")
  }

  # each sample get colname
  csv_dir <- output_dir
  csv_files <- dir(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  for (file in csv_files) {
    data <- read.csv(file, check.names = F)
    file_name <- tools::file_path_sans_ext(basename(file))  # 获取不带扩展名的文件名
    # 添加文件名前缀到列名（除了第一列）
    colnames(data)[-1] <- paste(file_name, colnames(data)[-1], sep = "")
    write.csv(data, file = file, row.names = FALSE)
    cat("Processed:", file_name, "\n")
  }

  # merge sample result
  file <- list.files(pattern = "*csv")
  data <- lapply(file, function(file) {
    df <- fread(file, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
    return(df)
  })
  data <- lapply(data, as.data.table)
  merged_data <- Reduce(function(x, y) merge(x, y, by = "gene_symbol", all = TRUE), data)
  df_selected <- merged_data %>%
    select(contains("_RE"))
  cols_to_select <- c("gene_symbol", colnames(df_selected))
  RE <- merged_data[, ..cols_to_select]
  setwd("../")
  write.csv(RE, file = paste0(out_dir,"3UTR_RE_Dealed.csv"), row.names = FALSE)


  ## IPA Quant ##
  dfIPA <- read.table(paste0(out_dir,"/updIPA_hg19.tsv"), header = TRUE)
  dfLE <- read.table(paste0(out_dir,"/updLE_hg19.tsv"), header = TRUE)

  # Quant each sample
  bamdir <- 'XXX'
  flsall <- dir(bamdir, pattern = ".Aligned.sortedByCoord.out.bam")
  flsall <- paste0(bamdir, flsall)  # 完整路径
  names(flsall) <- dir(bamdir, pattern = ".Aligned.sortedByCoord.out.bam")
  # 设置输出目录
  output_dir <- "/dell_1/yzh/apa/analysis/XXX/tcga/APAlyzer/Update_hg19/IPA/RE"

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  for (file in flsall) {
    identifier <- basename(file)
    IPA_OUTraw = PASEXP_IPA(dfIPA, dfLE, file, Strandtype = "forward", SeqType = "ThreeMostPairEnd", nts = 2)
    output <- file.path(output_dir, paste0(identifier, ".csv"))
    write.csv(IPA_OUTraw, file = output, row.names = FALSE, col.names = TRUE)
    cat("Processed:", identifier, "\n")
  }

  # each sample get colname
  csv_dir <- output_dir
  csv_files <- dir(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  for (file in csv_files) {
    data <- read.csv(file, check.names = F)
    file_name <- tools::file_path_sans_ext(basename(file))  # 获取不带扩展名的文件名
    # 添加文件名前缀到列名（除了前三列）
    colnames(data)[-c(1, 2, 3)] <- paste(file_name, colnames(data)[-c(1, 2, 3)], sep = "")
    write.csv(data, file = file, row.names = FALSE)
    cat("Processed:", file_name, "\n")
  }

  #merge sample result
  file <- list.files(pattern = "*csv")
  data <- lapply(file, function(file) {
    df <- fread(file, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
    return(df)
  })
  data <- lapply(data, as.data.table)
  merged_data <- Reduce(function(x, y) merge(x, y, by = c("IPAID", "gene_symbol", "PASid"), all = TRUE), data)
  df_selected <- merged_data %>%
    select(contains("_RE"))
  cols_to_select <- c("IPAID", "gene_symbol", "PASid", colnames(df_selected))
  RE <- merged_data[, ..cols_to_select]
  write.csv(RE, file = paste0(out_dir,"/IPA_RE_Dealed.csv"), row.names = FALSE)

}