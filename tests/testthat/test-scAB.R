library(Seurat)
library(dplyr)

# 创建测试数据生成函数
test_that("generate data", {
  # 1. 生成单细胞数据
  # 设置参数
  n_genes <- 100 # 100个基因
  n_cells <- 100 # 100个细胞
  n_cell_types <- 3 # 3种细胞类型

  # 生成基因名
  gene_names <- paste0("Gene", 1:n_genes)

  # 生成细胞类型标签
  cell_types <- paste0(
    "CellType",
    sample(1:n_cell_types, n_cells, replace = TRUE)
  )

  # 生成单细胞表达矩阵 (模拟UMI count数据)
  set.seed(42)
  sc_counts <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(gene_names, paste0("Cell", 1:n_cells))
  )

  # 为不同细胞类型添加一些差异表达信号
  for (i in 1:n_cell_types) {
    type_genes <- gene_names[(i - 1) * 30 + 1:30]
    type_cells <- which(cell_types == paste0("CellType", i))
    if (length(type_cells) > 0) {
      sc_counts[type_genes, type_cells] <- sc_counts[type_genes, type_cells] * 3
    }
  }

  # 2. 生成bulk数据 (基于单细胞数据聚合)
  n_bulk_samples <- 50 # 50个bulk样本

  # 创建bulk样本的细胞比例
  proportions <- matrix(
    runif(n_cell_types * n_bulk_samples),
    nrow = n_cell_types,
    ncol = n_bulk_samples
  )
  proportions <- proportions / rowSums(proportions) # 归一化

  # 初始化bulk表达矩阵
  bulk_counts <- matrix(
    0,
    nrow = n_genes,
    ncol = n_bulk_samples,
    dimnames = list(gene_names, paste0("Sample", 1:n_bulk_samples))
  )

  # 根据细胞比例生成bulk数据
  for (s in 1:n_bulk_samples) {
    for (i in 1:n_cell_types) {
      type_cells <- which(cell_types == paste0("CellType", i))
      if (length(type_cells) > 0) {
        # 获取该类型细胞的平均表达
        type_avg <- rowMeans(sc_counts[, type_cells])
        # 按照比例添加到bulk样本
        bulk_counts[, s] <- bulk_counts[, s] +
          type_avg * proportions[i, s] * 1000
      }
    }
    # 添加一些噪声
    bulk_counts[, s] <- rpois(n_genes, lambda = bulk_counts[, s])
  }

  # 3. 生成表型数据
  pheno_data <- data.frame(
    SampleID = paste0("Sample", 1:n_bulk_samples),
    Condition = sample(
      c("Control", "Treatment"),
      n_bulk_samples,
      replace = TRUE
    ),
    Age = sample(20:80, n_bulk_samples, replace = TRUE),
    stringsAsFactors = TRUE
  )

  # 确保表型数据与bulk样本匹配
  rownames(pheno_data) <- pheno_data$SampleID

  # 4. Seurat预处理
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = sc_counts, project = "TestData")

  # 添加细胞类型信息
  seurat_obj$cell_type <- cell_types

  # 标准Seurat预处理流程
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 50
  )
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(object = seurat_obj),
    npcs = 10
  )
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

  # 5. 验证数据处理结果
  # 验证Seurat对象处理正确
  expect_true("umap" %in% names(seurat_obj))
  expect_true(nrow(Embeddings(seurat_obj, reduction = "umap")) == n_cells)

  # 验证单细胞和bulk数据的基因匹配
  expect_true(all(rownames(sc_counts) %in% rownames(bulk_counts)))
  expect_true(
    length(intersect(rownames(sc_counts), rownames(bulk_counts))) == n_genes
  )

  # 验证表型数据与bulk样本匹配
  expect_true(all(rownames(pheno_data) %in% colnames(bulk_counts)))
  expect_true(
    length(intersect(rownames(pheno_data), colnames(bulk_counts))) ==
      n_bulk_samples
  )

  # 输出一些基本信息用于调试
  cat("\n=== Test Data Summary ===\n")
  cat(sprintf(
    "Single-cell data: %d genes, %d cells\n",
    nrow(sc_counts),
    ncol(sc_counts)
  ))
  cat(sprintf(
    "Bulk data: %d genes, %d samples\n",
    nrow(bulk_counts),
    ncol(bulk_counts)
  ))
  cat(sprintf("Phenotype data: %d samples\n", nrow(pheno_data)))
  cat("Cell types distribution:", table(cell_types), "\n")
  cat("Conditions distribution:", table(pheno_data$Condition), "\n")
  cat("=== End Summary ===\n")

  pheno <- setNames(
    ifelse(pheno_data$Condition == "Treatment", 1, 0),
    pheno_data$SampleID
  )

  # 保存处理后的数据以便进一步测试
  test_data <- list(
    single_cell = seurat_obj,
    bulk_data = bulk_counts,
    pheno_data = pheno
  )
})


test_that("scAB work", {
  scab_obj <- create_scAB.v5(
    Object = test_data$single_cell,
    bulk_dataset = test_data$bulk_data,
    phenotype = test_data$pheno_data,
    method = "binary",
    verbose = TRUE
  )

  expect_s3_class(scab_obj, "scAB_data")

  # This simulation data has negative values
  # But in real data, negative values are not allowed and indicate
  #   some technical issues.
  scab_obj$X[scab_obj$X < 0] <- 0

  k <- select_K.optimized(scab_obj, verbose = TRUE)

  expect_true(k > 0)

  loss_list <- scAB.optimized(Object = scab_obj, K = k)

  expect_equal(class(loss_list), "list")

  screened_seurat <- findSubset.optimized(
    Object = test_data$single_cell,
    scAB_Object = loss_list
  )

  expect_true("scAB" %in% colnames(screened_seurat[[]]))
})
