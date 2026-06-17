# ==============================================================================
# CondPED 批量模拟脚本
# 流程：生成数据 → 调用 CondPED 分析 → 保存指标（不保存原始数据）
# ==============================================================================

# --- 加载依赖 ---
library(CondPED)   # 你的 R 包（或 source R/ 目录下的文件）
source("inst/simulation/simulate_data.R")

# --- 配置 ---
n_rep <- 500
n <- 1000
p <- 1000

# 创建结果目录
dir.create("inst/simulation/results/sim1", recursive = TRUE, showWarnings = FALSE)
dir.create("inst/simulation/results/sim2", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Sim 1: 统计校准（参数扫描）
# ==============================================================================
run_sim1 <- function() {
  params <- expand.grid(
    n = c(200, 500, 1000, 2000),
    pve = c(0.01, 0.02, 0.05),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(params))) {
    n_i <- params$n[i]
    pve_i <- params$pve[i]
    cat(sprintf("\n=== Sim1: n=%d, PVE=%.2f%% ===\n", n_i, pve_i * 100))
    
    metrics <- vector("list", n_rep)
    
    for (r in 1:n_rep) {
      # 1. 生成数据
      dat <- generate_dataset_I(
        n = n_i, p = p,
        class1_pve = pve_i, class2_pve = pve_i,
        seed = r + i * 10000
      )
      
      # 2. 调用 CondPED 分析（假设函数已就绪）
      # result <- condped(X = dat$X, Y = dat$Y, alpha1 = NULL, alpha2 = 0.05)
      
      # 3. 提取指标（占位：先保存 truth 和方差验证）
      metrics[[r]] <- list(
        rep = r,
        n = n_i,
        pve = pve_i,
        var_y = apply(dat$Y, 2, var),           # 方差验证
        truth = dat$truth                         # 真值标签
        # 后续补充：
        # layer1_pvals = result$layer1$pvals,
        # detected = result$layer1$detected,
        # theta_cond = result$layer2$theta,
        # coverage = result$layer2$coverage
      )
      
      if (r %% 100 == 0) cat("  Rep", r, "done\n")
    }
    
    # 4. 保存该参数组合的结果
    save(metrics, file = sprintf(
      "inst/simulation/results/sim1/sim1_n%d_pve%d.RData",
      n_i, round(pve_i * 100)
    ))
  }
}

# ==============================================================================
# Sim 2: 核心分类验证（三个数据集各 500 次）
# ==============================================================================
run_sim2 <- function() {
  datasets <- list(
    I   = function(seed) generate_dataset_I(n = n, p = p, seed = seed),
    II  = function(seed) generate_dataset_II(n = n, p = p, seed = seed),
    III = function(seed) generate_dataset_III(n = n, p = p, seed = seed)
  )
  
  for (ds_name in names(datasets)) {
    cat(sprintf("\n=== Sim2 Dataset %s ===\n", ds_name))
    
    metrics <- vector("list", n_rep)
    seed_offset <- switch(ds_name, I = 0, II = 100000, III = 200000)
    
    for (r in 1:n_rep) {
      # 1. 生成数据
      dat <- datasets[[ds_name]](seed = r + seed_offset)
      
      # 2. 提取 IV 索引（Layer 3 需要）
      iv_idx <- which(dat$truth$is_IV)
      
      # 3. 调用 CondPED 分析（占位）
      # result <- condped(
      #   X = dat$X, Y = dat$Y,
      #   iv_idx = iv_idx,
      #   alpha1 = NULL,
      #   alpha2 = 0.05 / (p * 2)
      # )
      
      # 4. 保存指标（不保存原始 X/Y，太大）
      metrics[[r]] <- list(
        rep = r,
        dataset = ds_name,
        var_y = apply(dat$Y, 2, var),
        truth = dat$truth,
        Tau = dat$Tau
        # 后续补充：
        # predicted_class = result$predicted_class,
        # n_cond = result$layer2$n_cond,
        # tau_hat = result$layer3$tau_hat,
        # tau_pval = result$layer3$tau_pval
      )
      
      if (r %% 100 == 0) cat("  Rep", r, "done\n")
    }
    
    save(metrics, file = sprintf(
      "inst/simulation/results/sim2/sim2_dataset_%s.RData", ds_name
    ))
  }
}

# ==============================================================================
# 执行（后台跑时取消注释）
# ==============================================================================
# run_sim1()
# run_sim2()