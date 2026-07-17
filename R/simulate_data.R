# ==============================================================================
# simulate_data.R — CondPED 模拟数据生成模块
# ==============================================================================
#
# 【模块概述】
#   本文件实现 CondPED 论文中所有模拟数据集的生成函数。
#   共支持五类数据集：
#     Dataset I  ：无因果链（tau=0），验证 Class 1 / Class 2
#     Dataset II ：单向因果 A→B，验证 Class 1 / Class 3 / Class 4
#     Dataset III：双向因果 A↔B，验证 Class 5
#     Dataset IV ：三性状因果拓扑，验证 pairwise MR vs MVMR（补充模拟）
#     Dataset IV-misspec：四性状，真混杂 C + 诱饵 D（补充模拟）
#
# 【核心设计原则】
#   1. 功能位点采用固定效应（非随机），确保 Layer 1 稳定检出
#   2. IV位点数量 n_iv=20，PVE 动态设置为 max(0.02, 12/n)，保证 E[F]>12
#   3. 通过迭代调整残差协方差矩阵使 Var(Y_i)=1
#   4. 上位性和 G×E 互作默认关闭，在补充模拟 S5 中启用
#   5. Dataset II/III 中的 IV 位点在 truth 表里标记为 class="class3"：
#      生物学上 IV 和 Class 3 机制完全相同（完全中介），
#      is_IV=TRUE 保留，供 Layer 3 筛选工具变量时使用
#
# 【多环境扩展】
#   sim_phenotype_multienv() 支持 K 个环境，包含：
#     - 环境主效应 e_{ik} ~ N(0, sigma_env²)
#     - 基因型×环境互作 ae_{ilk} ~ N(0, sigma_gxe²)（随机斜率）
#   Layer 1 的 env 参数传入环境标识符后，自动切换为含环境固定效应的 OLS
#   Layer 2 的 method="LMM" 分支处理多环境条件效应估计
# ==============================================================================


# ==============================================================================
# 第一部分：基础工具函数
# ==============================================================================

#' 模拟基因型矩阵
#'
#' 根据指定群体类型生成 n×p 基因型矩阵，支持 RIL 和 F2 群体编码。
#'
#' 【编码方式】
#'   RIL (-1/+1)：纯合亲本1=-1，纯合亲本2=+1，由 maf 控制频率
#'   F2 (-1/0/+1)：纯合亲本1=-1，杂合子=0，纯合亲本2=+1
#'                 用二项分布(n=2, p=maf)生成后重编码
#'
#' 【基因型方差】
#'   RIL: Var(X) = 4 × maf × (1-maf)，maf=0.3时约为0.84
#'   F2:  Var(X) = 2 × maf × (1-maf)，maf=0.3时约为0.42
#'
#' @param n          整数。个体数（样本量）。
#' @param p          整数。SNP 数量。
#' @param maf        数值。次等位基因频率，默认0.3。
#' @param population 字符串。"RIL"（重组自交系）或"F2"（F2群体），默认"RIL"。
#' @param seed       整数或NULL。随机种子，保证结果可重复。
#' @return n×p 数值矩阵，行名为"Ind1...Indn"，列名为"SNP1...SNPp"。
#' @export
sim_genotype <- function(n, p,
                         maf = 0.3,
                         population = c("RIL", "F2"),
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  population <- match.arg(population)

  if (population == "RIL") {
    # RIL群体：每个位点为纯合子，以maf概率取值+1，以(1-maf)概率取值-1
    X <- matrix(
      sample(c(-1L, 1L), n * p, replace = TRUE, prob = c(1 - maf, maf)),
      nrow = n, ncol = p
    )
  } else {
    # F2群体：用二项分布模拟两个独立等位基因，然后重编码为-1/0/+1
    X <- matrix(rbinom(n * p, size = 2L, prob = maf), nrow = n, ncol = p)
    X[X == 0L] <- -1L # AA纯合子（次等位基因=0个）
    X[X == 1L] <- 0L # Aa杂合子（次等位基因=1个）
    X[X == 2L] <- 1L # aa纯合子（次等位基因=2个）
  }

  colnames(X) <- paste0("SNP", seq_len(p))
  rownames(X) <- paste0("Ind", seq_len(n))
  X
}


#' 由 PVE 计算加性效应量
#'
#' 在 Var(Y)=1 的标准化条件下，由位点解释方差比例（PVE）计算效应量 beta：
#'   beta = sqrt(PVE / Var(X))
#'
#' 这确保了无论群体类型如何，每个位点贡献的表型方差恰好等于设定的 PVE。
#'
#' @param pve        数值。该位点的表型方差解释比例（0到1之间）。
#' @param maf        数值。次等位基因频率，默认0.3。
#' @param population 字符串。"RIL"或"F2"，决定基因型方差的计算方式。
#' @return 数值标量。加性效应量 beta。
#' @export
beta_from_pve <- function(pve, maf = 0.3, population = c("RIL", "F2")) {
  population <- match.arg(population)
  # 根据群体类型计算基因型方差
  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  sqrt(pve / var_x)
}


# ==============================================================================
# 第二部分：方差标准化
# ==============================================================================

#' 计算使 Var(Y_i)=1 的残差协方差矩阵
#'
#' 【目的】确保所有性状的表型方差均为1，便于跨参数比较。
#'
#' 【数学原理】
#'   对于无因果链模型（M1，tau=0）：
#'     Sigma_Y = var_x × B'B + Sigma_E
#'     令 diag(Sigma_Y) = 1，可直接解析得 diag(Sigma_E) = 1 - diag(var_x×B'B)
#'
#'   对于有因果链模型（M2，tau≠0）：
#'     Sigma_Y = A' × (var_x×B'B + Sigma_E) × A
#'     其中 A = (I - tau')^{-1}
#'     需要迭代调整 diag(Sigma_E) 直到 diag(Sigma_Y) ≈ 1
#'
#' 【残差相关】rho>0时，性状间残差协方差 = rho × sqrt(sigma_i²×sigma_j²)
#'   这模拟共同环境因素（如温度、施肥）导致的跨性状表型相关
#'
#' @param B        p×m SNP效应矩阵。B[l,i]=位点l对性状i的加性效应。
#' @param Tau      m×m因果效应矩阵。Tau[j,i]=性状i对性状j的因果效应。
#' @param var_x    数值。基因型方差（RIL:4×maf×(1-maf)）。
#' @param rho      数值。残差相关系数，默认0（性状间残差独立）。
#' @param max_iter 整数。最大迭代次数，默认200。
#' @param tol      数值。收敛判断阈值，默认1e-8。
#' @return m×m 残差协方差矩阵 Sigma_E。
#' @export
compute_residual_variance <- function(B, Tau, var_x,
                                      rho = 0,
                                      max_iter = 200L,
                                      tol = 1e-8) {
  m <- ncol(B)
  Sigma_G <- var_x * crossprod(B) # 遗传协方差矩阵：var_x × B'B

  # ---- M1 闭合解（tau=0，无因果链）-----------------------------------------
  if (all(abs(Tau) < 1e-12)) {
    # 残差方差 = 1 - 遗传方差，保证Var(Y_i)=1
    diag_vals <- pmax(1e-6, 1 - diag(Sigma_G)) # 不低于1e-6防止奇异
    Sigma_E <- diag(diag_vals)

    # 若有残差相关，构造非对角元素
    if (abs(rho) > 1e-10 && m > 1L) {
      s <- sqrt(diag_vals)
      Sigma_E <- matrix(rho * outer(s, s), m, m)
      diag(Sigma_E) <- diag_vals
    }
    return(Sigma_E)
  }

  # ---- M2 迭代求解（tau≠0，有因果链）----------------------------------------
  # A = (I - Tau')^{-1}，是因果链传播算子
  A <- solve(diag(m) - t(Tau))

  # 初始估计：忽略因果链的影响
  Sigma_E <- diag(pmax(1e-6, 1 - diag(Sigma_G)))
  if (abs(rho) > 1e-10 && m > 1L) {
    s <- sqrt(diag(Sigma_E))
    Sigma_E <- matrix(rho * outer(s, s), m, m)
    diag(Sigma_E) <- s^2
  }

  # 迭代调整直到收敛
  for (iter in seq_len(max_iter)) {
    # 计算当前参数下的表型协方差
    Sigma_Y <- t(A) %*% (Sigma_G + Sigma_E) %*% A
    # 残差 = 目标方差1 - 当前方差
    diff <- 1 - diag(Sigma_Y)
    # 调整残差方差
    diag(Sigma_E) <- diag(Sigma_E) + diff

    # 更新残差相关（保持相关系数 rho 不变）
    if (abs(rho) > 1e-10 && m > 1L) {
      s <- sqrt(pmax(1e-6, diag(Sigma_E)))
      for (i in seq_len(m)) {
        for (j in seq_len(m)) {
          if (i != j) Sigma_E[i, j] <- rho * s[i] * s[j]
        }
      }
    }

    if (max(abs(diff)) < tol) break
    if (iter == max_iter) warning("方差标准化未收敛，请检查因果参数设置。")
  }
  Sigma_E
}


# ==============================================================================
# 第三部分：表型生成模型
# ==============================================================================

#' 单环境表型生成（模型M1，无因果链）
#'
#' 实现论文方程 M1：
#'   Y = X × B + eps,   eps ~ MVN(0, Sigma_E)
#'
#' 【上位性支持】
#'   若指定 epi_pairs 和 epi_effects，则在遗传主效应基础上叠加
#'   加性×加性上位性效应（aa项），对应论文 Eq.1 中的 aa_{ilh} 项。
#'   上位性编码：x_l × x_h（两个位点基因型的乘积）
#'
#' 【函数用途】
#'   本函数也被 sim_phenotype_M2 内部调用，生成"基础表型"后再传播因果链。
#'
#' @param X           n×p 基因型矩阵。
#' @param B           p×m 加性效应矩阵。B[l,i]=位点l对性状i的效应。
#' @param Sigma_E     m×m 残差协方差矩阵（由 compute_residual_variance 计算）。
#' @param epi_pairs   列表，每个元素为 c(l1, l2)，定义一对上位性位点。
#'   NULL=无上位性（默认）。
#' @param epi_effects 数值矩阵，行数=上位性对数，列数=m。
#'   元素[k,i]=第k对位点对性状i的上位性效应 aa_{l1l2,i}。
#' @param seed        整数或NULL。随机种子。
#' @return 命名列表：Y（n×m表型矩阵），G（n×m遗传值），eps（n×m残差）。
#' @export
sim_phenotype_M1 <- function(X, B, Sigma_E,
                             epi_pairs = NULL,
                             epi_effects = NULL,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  m <- ncol(B)

  # 计算加性遗传值：G = X × B（n×m矩阵）
  G <- X %*% B

  # 叠加上位性效应（若指定）
  if (!is.null(epi_pairs) && !is.null(epi_effects)) {
    for (k in seq_along(epi_pairs)) {
      pair <- epi_pairs[[k]]
      # 上位性项：x_l1 × x_l2 × aa_l1l2（列向量×行向量=外积）
      G <- G + (X[, pair[1L]] * X[, pair[2L]]) %*% epi_effects[k, , drop = FALSE]
    }
  }

  # 生成多变量正态残差
  eps <- MASS::mvrnorm(n, mu = rep(0, m), Sigma = Sigma_E)
  Y <- G + eps
  colnames(Y) <- paste0("Trait", seq_len(m))
  rownames(Y) <- rownames(X)
  list(Y = Y, G = G, eps = eps)
}


#' 单环境表型生成（模型M2，有因果链）
#'
#' 实现论文方程 M2：
#'   Y_base = X × B + eps
#'   Y      = Y_base × (I - Tau')^{-1}
#'
#' 【因果链传播机制】
#'   (I - Tau')^{-1} 是因果传播算子，将上游性状的遗传效应传递到下游性状。
#'   例如：A→B（tau_BA=0.3），B的表型 = 直接遗传效应 + 0.3 × A的表型
#'
#' 【稳态协方差】
#'   Sigma_Y = A' × Sigma_base × A，其中 A = (I - Tau')^{-1}
#'   compute_residual_variance 通过迭代确保 diag(Sigma_Y) = 1
#'
#' @param X           n×p 基因型矩阵。
#' @param B           p×m 加性效应矩阵。
#' @param Tau         m×m 因果效应矩阵。Tau[j,i]=性状i对性状j的因果效应。
#' @param Sigma_E     m×m 残差协方差（已由 compute_residual_variance 校准）。
#' @param epi_pairs   同 sim_phenotype_M1。
#' @param epi_effects 同 sim_phenotype_M1。
#' @param seed        整数或NULL。
#' @return 命名列表：Y，G，eps，Tau（原始），I_T_inv（因果传播算子）。
#' @export
sim_phenotype_M2 <- function(X, B, Tau, Sigma_E,
                             epi_pairs = NULL,
                             epi_effects = NULL,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- ncol(B)

  # 先生成无因果链的基础表型
  base <- sim_phenotype_M1(X, B, Sigma_E, epi_pairs, epi_effects, seed = NULL)

  # 因果链传播：Y = Y_base × (I - Tau')^{-1}
  I_T_inv <- solve(diag(m) - t(Tau))
  Y <- base$Y %*% I_T_inv

  rownames(Y) <- rownames(X)
  colnames(Y) <- colnames(base$Y)
  list(Y = Y, G = base$G, eps = base$eps, Tau = Tau, I_T_inv = I_T_inv)
}


# ==============================================================================
# 第四部分：多环境表型生成（补充模拟 S5）
# ==============================================================================

#' 多环境表型生成（含 G×E 互作）
#'
#' 实现论文补充方程中的多环境模型：
#'   y_{ijk} = mu_i + a_{il}×x_{jl}（加性主效，跨环境固定）
#'           + e_{ik}（环境主效应，性状i在环境k的随机偏移）
#'           + ae_{ilk}×x_{jl}（基因型×环境互作，随机斜率）
#'           + eps_{ijk}（残差）
#'
#' 【G×E 互作的含义】
#'   ae_{ilk} ~ N(0, sigma_gxe²)：同一位点在不同环境下效应不同。
#'   sigma_gxe=0 退化为单环境模型（纯加性，无互作）。
#'   sigma_gxe/sigma_env 的比值决定了 G×E 互作的相对重要性。
#'
#' 【输出格式】
#'   Y_long：(nK)×m 长格式，适合 Layer 1 LMM 扫描（lme4 输入）
#'   Y_by_env：n×m×K 三维数组，适合单环境分析或可视化
#'   Y_mean：n×m 跨环境均值，用于计算遗传力
#'
#' 【遗传力计算】
#'   h²_i = Var(G_i_main) / Var(Y_mean_i)
#'   这里用加性主效应方差除以跨环境均值的表型方差来近似遗传力。
#'
#' @param X          n×p 基因型矩阵（跨环境不变，RIL/F2群体特性）。
#' @param B          p×m 加性主效应矩阵（跨环境稳定的遗传效应）。
#' @param Sigma_E    m×m 残差协方差矩阵。
#' @param K          整数。环境数量，默认3。
#' @param sigma_env  数值。环境主效应标准差（e_{ik}的SD），默认0.3。
#'   值越大，环境间表型均值差异越大。
#' @param sigma_gxe  数值。G×E互作效应标准差（ae_{ilk}的SD），默认0.1。
#'   NULL或0表示无G×E互作（退化为可重复性模型）。
#' @param gxe_loci   整数向量。哪些位点有G×E互作效应（按列索引）。
#'   NULL=所有有非零B效应的位点（功能位点）都有G×E。
#' @param Tau        m×m 因果效应矩阵（NULL=无因果链，使用M1模型）。
#' @param seed       整数或NULL。
#' @return 命名列表：
#'   Y_long    (nK)×m，长格式表型矩阵，行名格式"IndX_EnvY"
#'   Y_by_env  n×m×K 三维数组
#'   Y_mean    n×m，跨环境均值
#'   X_long    (nK)×p，基因型矩阵（按环境重复K次）
#'   env_id    长度nK的环境标识符向量，格式"Env1","Env2"...
#'   ind_id    长度nK的个体标识符向量
#'   G_main    n×m，加性主效应遗传值（不含G×E）
#'   h2_obs    m维向量，各性状的近似遗传力
#'   K, n, sigma_env, sigma_gxe（记录实际使用的参数值）
#' @export
sim_phenotype_multienv <- function(X, B, Sigma_E,
                                   K = 3L,
                                   sigma_env = 0.3,
                                   sigma_gxe = 0.1,
                                   gxe_loci = NULL,
                                   Tau = NULL,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  m <- ncol(B)

  # ---- 确定有G×E互作的位点 --------------------------------------------------
  # 默认：所有有非零效应的位点都参与G×E互作（生物学上合理）
  if (is.null(gxe_loci)) {
    gxe_loci <- which(rowSums(B^2) > 1e-10)
  }

  # ---- 加性主效应（跨环境固定）------------------------------------------------
  # G_main = X × B，对所有环境相同
  G_main <- X %*% B # n×m

  # ---- 逐环境生成表型 ----------------------------------------------------------
  Y_by_env <- array(
    NA_real_,
    dim = c(n, m, K),
    dimnames = list(
      rownames(X),
      paste0("Trait", seq_len(m)),
      paste0("Env", seq_len(K))
    )
  )

  for (k in seq_len(K)) {
    # 环境主效应：e_{ik} ~ N(0, sigma_env²)
    # 对所有个体相同（环境整体偏移），但在性状间独立
    env_main <- rnorm(m, 0, sigma_env) # m维向量

    # G×E随机斜率：ae_{ilk} ~ N(0, sigma_gxe²)
    # 每个环境、每个有互作的位点、每个性状都独立抽样
    B_gxe <- matrix(0, p, m)
    if (!is.null(sigma_gxe) && sigma_gxe > 0 && length(gxe_loci) > 0) {
      B_gxe[gxe_loci, ] <- matrix(
        rnorm(length(gxe_loci) * m, 0, sigma_gxe),
        nrow = length(gxe_loci), ncol = m
      )
    }
    G_gxe <- X %*% B_gxe # n×m，本环境的G×E遗传效应

    # 残差：eps_{ijk} ~ MVN(0, Sigma_E)，性状间可能有相关
    eps <- MASS::mvrnorm(n, mu = rep(0, m), Sigma = Sigma_E)

    # 组合：表型 = 加性主效 + G×E效 + 残差，然后叠加环境偏移
    Y_k <- sweep(G_main + G_gxe + eps, 2, env_main, "+")

    # 若有因果链，通过因果传播算子传播
    if (!is.null(Tau)) {
      I_T_inv <- solve(diag(m) - t(Tau))
      Y_k <- Y_k %*% I_T_inv
    }

    colnames(Y_k) <- paste0("Trait", seq_len(m))
    rownames(Y_k) <- rownames(X)
    Y_by_env[, , k] <- Y_k
  }

  # ---- 转为长格式（nK×m）------------------------------------------------------
  # 将K个环境的表型矩阵按行堆叠
  # 行名格式："Ind1_Env1","Ind1_Env2",...,"InduN_EnvK"
  Y_long <- do.call(rbind, lapply(seq_len(K), function(k) Y_by_env[, , k]))
  env_id <- rep(paste0("Env", seq_len(K)), each = n) # 环境标识
  ind_id <- rep(rownames(X), times = K) # 个体标识
  X_long <- X[rep(seq_len(n), times = K), ] # 基因型按环境重复

  rownames(Y_long) <- paste0(ind_id, "_", env_id)
  rownames(X_long) <- paste0(ind_id, "_", env_id)

  # ---- 计算遗传力（基于跨环境均值）--------------------------------------------
  Y_mean <- apply(Y_by_env, c(1, 2), mean) # n×m，个体跨环境均值
  h2_obs <- apply(G_main, 2, var) / apply(Y_mean, 2, var)

  list(
    Y_long    = Y_long, # LMM扫描的主要输入
    Y_by_env  = Y_by_env, # 三维数组，便于可视化
    Y_mean    = Y_mean, # 跨环境均值，用于遗传力计算
    X_long    = X_long, # 对应的长格式基因型
    env_id    = env_id, # 环境标识符
    ind_id    = ind_id, # 个体标识符
    G_main    = G_main, # 加性主效应遗传值
    h2_obs    = h2_obs,
    K         = K,
    n         = n,
    sigma_env = sigma_env,
    sigma_gxe = sigma_gxe
  )
}


# ==============================================================================
# 第五部分：数据集生成函数（主模拟）
# ==============================================================================

#' Dataset I：无因果场景（tau=0），验证 Class 1 / Class 2
#'
#' 【数据集设计】
#'   生成模型：M1（无因果链）
#'     Y_A = X × beta_A + eps_A
#'     Y_B = X × beta_B + eps_B  （残差独立，rho=0）
#'
#' 【位点布局（p=1000）】
#'   SNP1-2  → Class 1（性状特异）：仅影响 Trait A，Trait B 无直接效应
#'   SNP3-4  → Class 2（水平多效）：独立影响 Trait A 和 Trait B
#'   SNP5+   → Null：无遗传效应
#'
#' 【效应量设计】
#'   同类内两个位点效应比例为1:0.8，避免联合OLS中的完全共线性
#'   两个位点同向（不取反），确保条件投影正常工作
#'
#' 【用途】
#'   Sim 1（Figure 1）的核心数据集，验证 Layer 1 的统计校准性质：
#'   Type I error（QQ图）、Power/FDR（曲线）、Coverage/RMSE（条件效应）
#'
#' @param n          整数。样本量。
#' @param p          整数。SNP数量，默认1000。
#' @param class1_pve 数值。Class 1位点对Trait A的PVE，默认0.05（5%）。
#' @param class2_pve 数值。Class 2位点对每个性状的PVE，默认0.05。
#' @param rho        数值。残差相关系数，默认0（两性状残差独立）。
#' @param epi_pairs  列表。NULL=无上位性。
#' @param epi_effects 矩阵。NULL=无上位性。
#' @param maf        数值。次等位基因频率，默认0.3。
#' @param population 字符串。"RIL"或"F2"，默认"RIL"。
#' @param seed       整数或NULL。
#' @return 命名列表：Y, G, eps, X, truth, dataset="I", seed, h2_obs。
#' @export
generate_dataset_I <- function(n, p = 1000L,
                               class1_pve = 0.05,
                               class2_pve = 0.05,
                               rho = 0,
                               epi_pairs = NULL,
                               epi_effects = NULL,
                               maf = 0.3,
                               population = "RIL",
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- sim_genotype(n, p, maf, population)
  B <- matrix(0, p, 2L)
  b1 <- beta_from_pve(class1_pve, maf, population)
  b2 <- beta_from_pve(class2_pve, maf, population)

  # Class 1：SNP1-2，仅影响Trait A，同向效应，比例1:0.8
  B[1L, 1L] <- b1
  B[2L, 1L] <- b1 * 0.8

  # Class 2：SNP3-4，同时独立影响Trait A和B，各自同向
  B[3L, 1L] <- b2
  B[3L, 2L] <- b2
  B[4L, 1L] <- b2 * 0.8
  B[4L, 2L] <- b2 * 0.8

  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, matrix(0, 2L, 2L), var_x, rho = rho)
  out <- sim_phenotype_M1(X, B, Sigma_E, epi_pairs, epi_effects)

  if (is.null(rownames(out$Y))) rownames(out$Y) <- rownames(X)
  if (is.null(colnames(out$Y))) colnames(out$Y) <- paste0("Trait", seq_len(2L))
  stopifnot(identical(rownames(out$Y), rownames(X)))

  h2_obs <- diag(var(out$G)) / diag(var(out$Y))

  truth <- data.frame(
    SNP = colnames(X),
    class = c(rep("class1", 2L), rep("class2", 2L), rep("null", p - 4L)),
    beta_A = B[, 1L],
    beta_B = B[, 2L],
    is_IV = FALSE,
    stringsAsFactors = FALSE
  )

  c(out, list(X = X, truth = truth, dataset = "I", seed = seed, h2_obs = h2_obs))
}


#' Dataset II：单向因果 A→B，验证 Class 1 / Class 3 / Class 4
#'
#' 【数据集设计】
#'   生成模型：M2（单向因果链，A→B）
#'     Y_A = X × beta_A + eps_A
#'     Y_B = tau × Y_A + X × beta_B_direct + eps_B
#'
#' 【位点布局】
#'   SNP1-2      → Class 1：仅影响Trait B（beta_A=0）
#'   SNP3-4      → Class 3（完全中介）：仅影响Trait A，B效应完全经因果链传递
#'   SNP5-6      → Class 4（部分中介）：同时影响A（强效应）和B（弱直接效应）
#'   SNP7-SNP(6+n_iv) → IV辅助位点：仅影响A，用于Layer 3 MR估计tau
#'   SNP(7+n_iv)+ → Null
#'
#' 【IV位点标签说明（关键设计决策）】
#'   IV位点生物学机制：仅影响上游性状A，B的关联完全通过因果链A→B传递
#'   这与Class 3完全相同，因此truth表中class="class3"
#'   is_IV=TRUE保留，供select_instruments()筛选工具变量时识别
#'   效应量略小（PVE=max(2%,12/n)）以与主要功能位点区分，但机制相同
#'
#' 【IV PVE动态设置】
#'   iv_pve = max(0.02, 12/n)，确保期望F统计量 E[F] > 12
#'   F = n × PVE / (1-PVE)，F>10是弱工具变量的经典判断标准
#'   n=200时：max(0.02, 0.06)=0.06；n=1000时：max(0.02, 0.012)=0.02
#'
#' @param n              整数。样本量。
#' @param p              整数。SNP数量，默认1000。
#' @param class1_pve     数值。Class 1位点对Trait B的PVE，默认0.05。
#' @param class3_pve     数值。Class 3位点对Trait A的PVE，默认0.05。
#' @param class4_pve_A   数值。Class 4位点对Trait A的PVE，默认0.05。
#' @param class4_pve_B   数值。Class 4位点对Trait B的直接效应PVE，默认0.05。
#'   降低此值（如0.01）可测试弱部分中介的检出能力。
#' @param tau            数值。A→B的因果效应，默认0.3。
#' @param iv_pve         数值或NULL。IV位点PVE。NULL=自动计算max(0.02,12/n)。
#' @param n_iv           整数。IV辅助位点数量，默认20。
#' @param rho            数值。残差相关系数，默认0。
#' @param epi_pairs,epi_effects  NULL=无上位性。
#' @param maf,population  同其他函数。
#' @param seed           整数或NULL。
#' @return 命名列表：Y, G, eps, X, truth, dataset="II", Tau, seed, h2_obs, iv_pve_used。
#' @export
generate_dataset_II <- function(n, p = 1000L,
                                class1_pve = 0.05,
                                class3_pve = 0.05,
                                class4_pve_A = 0.05,
                                class4_pve_B = 0.05,
                                tau = 0.3,
                                iv_pve = NULL,
                                n_iv = 20L,
                                rho = 0,
                                epi_pairs = NULL,
                                epi_effects = NULL,
                                maf = 0.3,
                                population = "RIL",
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 动态计算IV PVE，保证工具变量强度（E[F]>12）
  if (is.null(iv_pve)) iv_pve <- max(0.02, 12 / n)

  n_iv <- as.integer(n_iv)
  n_func <- 6L + n_iv # 总功能位点数（不含Null）
  if (n_func >= p) stop(sprintf("n_iv=%d 过大（p=%d）。", n_iv, p))

  X <- sim_genotype(n, p, maf, population)
  B <- matrix(0, p, 2L)

  # 计算各类别的效应量
  b1 <- beta_from_pve(class1_pve, maf, population) # Class 1对B
  b3 <- beta_from_pve(class3_pve, maf, population) # Class 3对A
  b4a <- beta_from_pve(class4_pve_A, maf, population) # Class 4对A
  b4b <- beta_from_pve(class4_pve_B, maf, population) # Class 4对B（直接）
  biv <- beta_from_pve(iv_pve, maf, population) # IV对A

  # 设置各类位点的效应：同向，比例1:0.8
  B[1L, 2L] <- b1
  B[2L, 2L] <- b1 * 0.8 # Class 1：B特异性
  B[3L, 1L] <- b3
  B[4L, 1L] <- b3 * 0.8 # Class 3：A特异性（B完全中介）
  B[5L, 1L] <- b4a
  B[5L, 2L] <- b4b # Class 4：A强效+B弱直接效应
  B[6L, 1L] <- b4a * 0.8
  B[6L, 2L] <- b4b * 0.8

  # IV位点：仅影响A，交替正负（增加多样性，避免LD修剪时全部被滤除）
  iv_idx <- 7L:(6L + n_iv)
  B[iv_idx, 1L] <- biv * rep(c(1, -1), length.out = n_iv)

  # 设置因果矩阵：A→B（tau_BA=tau）
  Tau <- matrix(0, 2L, 2L)
  Tau[2L, 1L] <- tau

  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, Tau, var_x, rho = rho)
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E, epi_pairs, epi_effects)

  if (is.null(rownames(out$Y))) rownames(out$Y) <- rownames(X)
  if (is.null(colnames(out$Y))) colnames(out$Y) <- paste0("Trait", seq_len(2L))
  stopifnot(identical(rownames(out$Y), rownames(X)))

  h2_obs <- diag(var(out$G)) / diag(var(out$Y))

  truth <- data.frame(
    SNP = colnames(X),
    # ── 关键设计说明 ───────────────────────────────────────────────────────
    # IV位点与Class 3在统计指纹上完全相同：
    #   仅影响Trait A，Trait B的关联完全由因果链A→B传递（beta_B=0）
    # 因此truth$class统一标记为"class3"，消除评估时的虚假MISMATCH
    # is_IV=TRUE保留，Layer 3中select_instruments()通过此列识别工具变量候选
    # ─────────────────────────────────────────────────────────────────────
    class = c(
      rep("class1", 2L),
      rep("class3", 2L),
      rep("class4", 2L),
      rep("class3", n_iv), # IV位点标记为class3（不再是"IV_A"）
      rep("null", p - 6L - n_iv)
    ),
    beta_A = B[, 1L],
    beta_B = B[, 2L],
    # tau列：标记哪些位点与A→B因果效应相关（IV位点作为该效应的工具变量）
    tau = c(rep(0, 6L), rep(tau, n_iv), rep(0, p - 6L - n_iv)),
    # is_IV列：标记工具变量位点，供MR筛选使用
    is_IV = c(rep(FALSE, 6L), rep(TRUE, n_iv), rep(FALSE, p - 6L - n_iv)),
    stringsAsFactors = FALSE
  )

  c(out, list(
    X = X, truth = truth, dataset = "II", Tau = Tau,
    seed = seed, h2_obs = h2_obs, iv_pve_used = iv_pve
  ))
}


#' Dataset III：双向因果 A↔B，验证 Class 5
#'
#' 【数据集设计】
#'   生成模型：M2（双向因果链，联立方程组）
#'     Y_A = tau_BA × Y_B + X × beta_A + eps_A
#'     Y_B = tau_AB × Y_A + X × beta_B + eps_B
#'
#' 【位点布局】
#'   SNP1-2           → Class 5（双向因果）：同时对A和B有直接效应
#'   SNP3-SNP(2+n_iv) → IV_A（A方向工具变量）：仅影响A，用于估计tau_AB
#'   SNP(3+n_iv)-SNP(2+2n_iv) → IV_B（B方向工具变量）：仅影响B，用于估计tau_BA
#'   SNP(3+2n_iv)+    → Null
#'
#' 【IV位点标签说明】
#'   IV_A位点：仅影响A，B的关联通过A→B传递 → 本质上是A→B方向的Class 3
#'   IV_B位点：仅影响B，A的关联通过B→A传递 → 本质上是B→A方向的Class 3
#'   两者均标记为class="class3"，is_IV=TRUE保留供MR使用
#'
#' @param n,p         样本量和SNP数。
#' @param class5_pve  数值。Class 5位点每个性状的PVE，默认0.05。
#' @param tau_AB      数值。A→B的因果效应，默认0.2。
#' @param tau_BA      数值。B→A的因果效应，默认0.2。
#' @param iv_pve      数值或NULL。IV位点PVE，NULL=自动。
#' @param n_iv        整数。每个方向的IV位点数，默认20。
#' @param rho,epi_pairs,epi_effects,maf,population,seed  同其他函数。
#' @return 命名列表：Y, G, eps, X, truth, dataset="III", Tau, seed, h2_obs, iv_pve_used。
#' @export
generate_dataset_III <- function(n, p = 1000L,
                                 class5_pve = 0.05,
                                 tau_AB = 0.2,
                                 tau_BA = 0.2,
                                 iv_pve = NULL,
                                 n_iv = 20L,
                                 rho = 0,
                                 epi_pairs = NULL,
                                 epi_effects = NULL,
                                 maf = 0.3,
                                 population = "RIL",
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(iv_pve)) iv_pve <- max(0.02, 12 / n)

  n_iv <- as.integer(n_iv)
  n_func <- 2L + 2L * n_iv
  if (n_func >= p) stop(sprintf("2×n_iv=%d 过大（p=%d）。", 2L * n_iv, p))

  X <- sim_genotype(n, p, maf, population)
  B <- matrix(0, p, 2L)
  b5 <- beta_from_pve(class5_pve, maf, population)
  biv <- beta_from_pve(iv_pve, maf, population)

  # Class 5：SNP1-2，同时影响A和B，同向，比例1:0.8
  B[1L, 1L] <- b5
  B[1L, 2L] <- b5
  B[2L, 1L] <- b5 * 0.8
  B[2L, 2L] <- b5 * 0.8

  # IV_A：仅影响A，用于估计A→B因果效应
  iv_A_idx <- 3L:(2L + n_iv)
  B[iv_A_idx, 1L] <- biv * rep(c(1, -1), length.out = n_iv)

  # IV_B：仅影响B，用于估计B→A因果效应
  iv_B_idx <- (3L + n_iv):(2L + 2L * n_iv)
  B[iv_B_idx, 2L] <- biv * rep(c(1, -1), length.out = n_iv)

  # 设置双向因果矩阵
  Tau <- matrix(0, 2L, 2L)
  Tau[2L, 1L] <- tau_AB # A → B
  Tau[1L, 2L] <- tau_BA # B → A

  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, Tau, var_x, rho = rho)
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E, epi_pairs, epi_effects)

  if (is.null(rownames(out$Y))) rownames(out$Y) <- rownames(X)
  if (is.null(colnames(out$Y))) colnames(out$Y) <- paste0("Trait", seq_len(2L))
  stopifnot(identical(rownames(out$Y), rownames(X)))

  h2_obs <- diag(var(out$G)) / diag(var(out$Y))

  truth <- data.frame(
    SNP = colnames(X),
    # IV_A：A方向完全中介，机制等同Class 3（A→B方向）
    # IV_B：B方向完全中介，机制等同Class 3（B→A方向）
    # 均标记为class="class3"，is_IV=TRUE区分二者
    class = c(
      rep("class5", 2L),
      rep("class3", n_iv), # IV_A → class3（原"IV_A"）
      rep("class3", n_iv), # IV_B → class3（原"IV_B"）
      rep("null", p - 2L - 2L * n_iv)
    ),
    beta_A = B[, 1L],
    beta_B = B[, 2L],
    # tau_AB列：IV_A位点是A→B方向的工具变量
    tau_AB = c(
      rep(tau_AB, 2L + n_iv), rep(0, n_iv), rep(0, p - 2L - 2L * n_iv)
    ),
    # tau_BA列：IV_B位点是B→A方向的工具变量
    tau_BA = c(
      rep(tau_BA, 2L), rep(0, n_iv), rep(tau_BA, n_iv), rep(0, p - 2L - 2L * n_iv)
    ),
    is_IV = c(
      rep(FALSE, 2L),
      rep(TRUE, n_iv), # IV_A
      rep(TRUE, n_iv), # IV_B
      rep(FALSE, p - 2L - 2L * n_iv)
    ),
    stringsAsFactors = FALSE
  )

  c(out, list(
    X = X, truth = truth, dataset = "III", Tau = Tau,
    seed = seed, h2_obs = h2_obs, iv_pve_used = iv_pve
  ))
}


# ==============================================================================
# 第六部分：数据集生成函数（补充模拟）
# ==============================================================================

#' Dataset IV：三性状因果拓扑，验证 pairwise MR vs MVMR（补充模拟）
#'
#' 【验证目的】
#'   对比三种因果拓扑下，两两MR和MVMR的估计偏差：
#'   1. chain（因果链 A→B→C）：无共同上游，两两MR理论无偏，验证良性情形
#'   2. confound_obs（C→A，C→B，真值tau_AB=0）：可观测混杂，两两MR产生假因果
#'   3. confound_latent（U→A，U→B，U不可测）：潜在混杂，两种方法均无法纠正
#'
#' 【位点布局（按拓扑不同）】
#'   chain：
#'     IV_A（n_iv个）：仅影响A（清洁工具变量）
#'     IV_B（n_iv个）：仅影响B（清洁工具变量）
#'     Cspecific（5个）：仅影响C
#'   confound_obs：
#'     IV_C_confounder（n_conf个）：直接影响混杂因子C（污染工具变量）
#'     IV_A_clean（n_iv个）：清洁工具变量
#'     IV_B_clean（n_iv个）：清洁工具变量
#'   confound_latent：
#'     IV_U_latent（n_conf个）：影响潜在混杂U（不可观测）
#'     IV_A_clean（n_iv个）：清洁工具变量
#'     IV_B_clean（n_iv个）：清洁工具变量
#'     Cspecific（5个）：仅影响C
#'
#' 【安全约束】
#'   n_conf × conf_pve < 1，否则方差标准化失败（残差方差为负）
#'
#' @param n          整数。样本量（保持与主模拟一致）。
#' @param p          整数。SNP数，默认1000。
#' @param topology   字符串。"chain"/"confound_obs"/"confound_latent"。
#' @param tau_AB     数值。chain拓扑中A→B因果效应，默认0.4。
#' @param tau_BC     数值。chain拓扑中B→C因果效应，默认0.4。
#' @param tau_CA     数值。confound拓扑中C→A的效应，默认0.5。
#' @param tau_CB     数值。confound拓扑中C→B的效应，默认0.5。
#' @param tau_UA,tau_UB  数值。latent拓扑中U→A和U→B的效应，默认0.5。
#' @param conf_pve   数值。混杂工具变量单点PVE，默认0.03。
#'   n_conf × conf_pve必须小于0.9。
#' @param class_pve  数值。性状特异位点PVE，默认0.05。
#' @param iv_pve     数值或NULL。清洁工具变量PVE，NULL=自动。
#' @param n_iv       整数。清洁工具变量数，默认20。
#' @param n_conf     整数。混杂工具变量数，默认10。
#' @param rho,maf,population,seed  同其他函数。
#' @return 命名列表：Y(n×3), G, eps, X, truth, dataset="IV", topology,
#'   Tau_full, Tau_obs(3×3), tau_AB_true, latent(逻辑值), seed, h2_obs, iv_pve_used。
#' @export
generate_dataset_IV <- function(n,
                                p = 1000L,
                                topology = c("chain", "confound_obs", "confound_latent"),
                                tau_AB = 0.4,
                                tau_BC = 0.4,
                                tau_CA = 0.5,
                                tau_CB = 0.5,
                                tau_UA = 0.5,
                                tau_UB = 0.5,
                                conf_pve = 0.03,
                                class_pve = 0.05,
                                iv_pve = NULL,
                                n_iv = 20L,
                                n_conf = 10L,
                                rho = 0,
                                maf = 0.3,
                                population = "RIL",
                                seed = NULL) {
  topology <- match.arg(topology)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(iv_pve)) iv_pve <- max(0.02, 12 / n)
  n_iv <- as.integer(n_iv)
  n_conf <- as.integer(n_conf)

  # 安全检查：防止方差标准化失败
  if (n_conf * conf_pve >= 0.9) {
    warning(sprintf(
      "n_conf×conf_pve=%.2f 接近1，方差标准化可能失败，建议降低n_conf或conf_pve。",
      n_conf * conf_pve
    ))
  }

  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  biv <- beta_from_pve(iv_pve, maf, population) # 清洁工具变量效应量
  bcf <- beta_from_pve(conf_pve, maf, population) # 混杂工具变量效应量
  bcl <- beta_from_pve(class_pve, maf, population) # 性状特异位点效应量
  sgn <- function(k) rep(c(1, -1), length.out = k) # 交替正负编码

  role <- rep("null", p) # 位点角色标签
  iv_for <- rep(NA_character_, p) # 工具变量指向的性状

  # ---- 按拓扑构建效应矩阵和因果矩阵 ------------------------------------------
  if (topology == "chain") {
    # 因果链 A→B→C：每个有向对之间有单向因果，无共同上游混杂
    m_full <- 3L
    n_func <- 2L * n_iv + 5L
    if (n_func >= p) stop(sprintf("位点块总数(%d)超过p(%d)。", n_func, p))

    B_full <- matrix(0, p, m_full)
    Tau_full <- matrix(0, m_full, m_full)
    Tau_full[2L, 1L] <- tau_AB # A → B
    Tau_full[3L, 2L] <- tau_BC # B → C

    # 布置位点
    iA <- 1L:n_iv
    B_full[iA, 1L] <- biv * sgn(n_iv) # A的工具变量
    iB <- (n_iv + 1L):(2L * n_iv)
    B_full[iB, 2L] <- biv * sgn(n_iv) # B的工具变量
    iC <- (2L * n_iv + 1L):(2L * n_iv + 5L)
    B_full[iC, 3L] <- bcl # C特异位点

    role[iA] <- "IV_A"
    iv_for[iA] <- "A"
    role[iB] <- "IV_B"
    iv_for[iB] <- "B"
    role[iC] <- "Cspecific"
  } else if (topology == "confound_obs") {
    # 可观测混杂：C→A，C→B，真值tau_AB=0
    # 作用于C的工具变量会经C→A污染A的工具集，再经C→B影响B，违反排他性
    m_full <- 3L
    n_func <- n_conf + 2L * n_iv
    if (n_func >= p) stop(sprintf("位点块总数(%d)超过p(%d)。", n_func, p))

    B_full <- matrix(0, p, m_full)
    Tau_full <- matrix(0, m_full, m_full)
    Tau_full[1L, 3L] <- tau_CA # C → A（混杂路径）
    Tau_full[2L, 3L] <- tau_CB # C → B（混杂路径）
    # Tau_full[2L,1L]=0：A→B真实效应为0（待估计的是假因果）

    iC <- 1L:n_conf
    B_full[iC, 3L] <- bcf * sgn(n_conf) # 混杂工具变量
    iA <- (n_conf + 1L):(n_conf + n_iv)
    B_full[iA, 1L] <- biv * sgn(n_iv) # 清洁工具变量
    iB <- (n_conf + n_iv + 1L):(n_conf + 2L * n_iv)
    B_full[iB, 2L] <- biv * sgn(n_iv)

    role[iC] <- "IV_C_confounder"
    iv_for[iC] <- "C"
    role[iA] <- "IV_A_clean"
    iv_for[iA] <- "A"
    role[iB] <- "IV_B_clean"
    iv_for[iB] <- "B"
  } else {
    # 潜在混杂：U→A，U→B，U不可测（m_full=4，但Y只保留前3列）
    # 两两MR和MVMR均无法纠正，表现为双向假信号（Class 5误判）
    m_full <- 4L
    n_func <- n_conf + 2L * n_iv + 5L
    if (n_func >= p) stop(sprintf("位点块总数(%d)超过p(%d)。", n_func, p))

    B_full <- matrix(0, p, m_full)
    Tau_full <- matrix(0, m_full, m_full)
    Tau_full[1L, 4L] <- tau_UA # U → A（第4列为潜在变量U）
    Tau_full[2L, 4L] <- tau_UB # U → B

    iU <- 1L:n_conf
    B_full[iU, 4L] <- bcf * sgn(n_conf) # 作用于潜在U
    iA <- (n_conf + 1L):(n_conf + n_iv)
    B_full[iA, 1L] <- biv * sgn(n_iv)
    iB <- (n_conf + n_iv + 1L):(n_conf + 2L * n_iv)
    B_full[iB, 2L] <- biv * sgn(n_iv)
    iC <- (n_conf + 2L * n_iv + 1L):(n_conf + 2L * n_iv + 5L)
    B_full[iC, 3L] <- bcl

    role[iU] <- "IV_U_latent"
    iv_for[iU] <- "U"
    role[iA] <- "IV_A_clean"
    iv_for[iA] <- "A"
    role[iB] <- "IV_B_clean"
    iv_for[iB] <- "B"
    role[iC] <- "Cspecific"
  }

  # ---- 生成表型（在完整m_full性状系统上），然后只保留前3个可观测性状 ----------
  X <- sim_genotype(n, p, maf, population)
  Sigma_E <- compute_residual_variance(B_full, Tau_full, var_x, rho = rho)
  out <- sim_phenotype_M2(X, B_full, Tau_full, Sigma_E)

  obs <- 1L:3L # 只保留前3个可观测性状（潜在变量U对应第4列，不输出）
  Y <- out$Y[, obs, drop = FALSE]
  G <- out$G[, obs, drop = FALSE]
  eps <- out$eps[, obs, drop = FALSE]
  colnames(Y) <- c("TraitA", "TraitB", "TraitC")
  rownames(Y) <- rownames(X)
  stopifnot(identical(rownames(Y), rownames(X)))

  h2_obs <- diag(var(G)) / diag(var(Y))

  truth <- data.frame(
    SNP = colnames(X),
    role = role,
    iv_for = iv_for,
    beta_A = B_full[, 1L],
    beta_B = B_full[, 2L],
    beta_C = B_full[, 3L],
    class = role, # 兼容check_dataset()
    is_IV = !is.na(iv_for),
    stringsAsFactors = FALSE
  )
  if (m_full == 4L) truth$beta_U <- B_full[, 4L] # 潜在变量效应（latent拓扑）

  Tau_obs <- Tau_full[obs, obs, drop = FALSE] # 可观测性状间的因果矩阵

  c(
    list(Y = Y, G = G, eps = eps),
    list(
      X = X, truth = truth, dataset = "IV", topology = topology,
      Tau_full = Tau_full,
      Tau_obs = Tau_obs,
      tau_AB_true = Tau_full[2L, 1L], # confound拓扑中为0，chain拓扑中为tau_AB
      latent = (m_full == 4L), # 是否有潜在变量
      seed = seed, h2_obs = h2_obs, iv_pve_used = iv_pve
    )
  )
}


#' Dataset IV-misspec：四性状，验证MVMR错误指定的失效（补充模拟）
#'
#' 【验证目的】
#'   证明MVMR的纠偏能力完全依赖正确指定混杂因子：
#'   - MVMR正确指定C（真混杂）→ Type I error ≈ 5%（无偏）
#'   - MVMR错误指定D（诱饵，是C的下游）→ Type I error ≈ 97%（与两两MR相当）
#'
#' 【性状角色】
#'   A=1：暴露性状（待估计A→B）
#'   B=2：结局性状（tau_AB=0，无真实因果）
#'   C=3：真混杂（C→A，C→B，产生A和B的假相关）
#'   D=4：诱饵（C→D，和系统相关但不是A/B的共同上游）
#'
#' 【关键参数安全约束】n_conf × conf_pve < 0.9
#'
#' @param n,p           样本量和SNP数。
#' @param tau_CA        数值。C→A效应，默认0.5。
#' @param tau_CB        数值。C→B效应，默认0.5。
#' @param tau_CD        数值。C→D效应（诱饵），默认0.5。
#' @param conf_pve      数值。混杂工具变量单点PVE，默认0.04。
#' @param iv_pve        数值或NULL。清洁工具变量PVE，NULL=自动。
#' @param n_iv          整数。每个清洁块的工具变量数，默认20。
#' @param n_conf        整数。混杂工具变量数，默认15。
#' @param rho,maf,population,seed  同其他函数。
#' @return 命名列表：Y(n×4), G, eps, X, truth, dataset="IV_misspec", Tau_full,
#'   tau_AB_true=0, seed, h2_obs, iv_pve_used。
#' @export
generate_dataset_IV_misspec <- function(n, p = 1000L,
                                        tau_CA = 0.5,
                                        tau_CB = 0.5,
                                        tau_CD = 0.5,
                                        conf_pve = 0.04,
                                        iv_pve = NULL,
                                        n_iv = 20L,
                                        n_conf = 15L,
                                        rho = 0,
                                        maf = 0.3,
                                        population = "RIL",
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(iv_pve)) iv_pve <- max(0.02, 12 / n)
  n_iv <- as.integer(n_iv)
  n_conf <- as.integer(n_conf)
  if (n_conf * conf_pve >= 0.9) {
    warning(sprintf("n_conf×conf_pve=%.2f 接近1，方差标准化可能失败。", n_conf * conf_pve))
  }

  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  biv <- beta_from_pve(iv_pve, maf, population)
  bcf <- beta_from_pve(conf_pve, maf, population)
  sgn <- function(k) rep(c(1, -1), length.out = k)

  m_full <- 4L
  n_func <- n_conf + 3L * n_iv
  if (n_func >= p) stop(sprintf("位点块总数(%d)超过p(%d)。", n_func, p))

  B <- matrix(0, p, m_full)
  Tau <- matrix(0, m_full, m_full)
  Tau[1L, 3L] <- tau_CA # C → A（真混杂路径）
  Tau[2L, 3L] <- tau_CB # C → B（真混杂路径）
  Tau[4L, 3L] <- tau_CD # C → D（诱饵：D是C的下游，非A/B的共同上游）
  # Tau[2L,1L]=0：A→B真实效应为0（所有A→B关联均为混杂导致的假因果）

  role <- rep("null", p)
  iv_for <- rep(NA_character_, p)

  # 位点布局：混杂工具→A清洁工具→B清洁工具→D清洁工具
  iC <- 1L:n_conf
  B[iC, 3L] <- bcf * sgn(n_conf) # 直接影响C
  iA <- (n_conf + 1L):(n_conf + n_iv)
  B[iA, 1L] <- biv * sgn(n_iv) # 清洁A工具
  iB <- (n_conf + n_iv + 1L):(n_conf + 2L * n_iv)
  B[iB, 2L] <- biv * sgn(n_iv) # 清洁B工具
  iD <- (n_conf + 2L * n_iv + 1L):(n_conf + 3L * n_iv)
  B[iD, 4L] <- biv * sgn(n_iv) # D工具（诱饵）

  role[iC] <- "IV_C_confounder"
  iv_for[iC] <- "C"
  role[iA] <- "IV_A_clean"
  iv_for[iA] <- "A"
  role[iB] <- "IV_B_clean"
  iv_for[iB] <- "B"
  role[iD] <- "IV_D_decoy"
  iv_for[iD] <- "D"

  X <- sim_genotype(n, p, maf, population)
  Sigma_E <- compute_residual_variance(B, Tau, var_x, rho = rho)
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E)

  # 四个性状全部可观测（与Dataset IV不同，此处无潜在变量）
  Y <- out$Y
  colnames(Y) <- c("TraitA", "TraitB", "TraitC", "TraitD")
  rownames(Y) <- rownames(X)
  h2_obs <- diag(var(out$G)) / diag(var(Y))

  truth <- data.frame(
    SNP = colnames(X), role = role, iv_for = iv_for,
    beta_A = B[, 1L], beta_B = B[, 2L], beta_C = B[, 3L], beta_D = B[, 4L],
    class = role, is_IV = !is.na(iv_for),
    stringsAsFactors = FALSE
  )

  list(
    Y = Y, G = out$G, eps = out$eps, X = X, truth = truth,
    dataset = "IV_misspec",
    Tau_full = Tau,
    Tau_obs = Tau,
    tau_AB_true = Tau[2L, 1L], # = 0，验证FP率
    seed = seed,
    h2_obs = h2_obs,
    iv_pve_used = iv_pve
  )
}


# ==============================================================================
# 第七部分：统一入口和诊断工具
# ==============================================================================

#' 统一数据集生成入口（支持 Dataset I / II / III / IV）
#'
#' 根据 dataset 参数调用对应的生成函数，所有额外参数通过 ... 透传。
#'
#' @param dataset 字符串。"I"/"II"/"III"/"IV"之一。
#'   Dataset IV 需额外传入 topology 参数（见 generate_dataset_IV）。
#' @param ...     传递给对应数据集生成函数的额外参数。
#' @return 对应生成函数的命名列表输出。
#' @export
generate_condped <- function(dataset = c("I", "II", "III", "IV"), ...) {
  dataset <- match.arg(dataset)
  switch(dataset,
    I   = generate_dataset_I(...),
    II  = generate_dataset_II(...),
    III = generate_dataset_III(...),
    IV  = generate_dataset_IV(...)
  )
}


#' 数据集快速诊断（检查生成结果是否正常）
#'
#' 打印表型方差、遗传力、位点计数和IV强度，帮助验证参数设置。
#' 期望结果：Trait variances ≈ 1，IV E[F] > 12。
#'
#' 【注意】
#'   功能位点计数：class∈{class1,...,class5}且is_IV=FALSE
#'   IV位点：is_IV=TRUE（包含标记为class3的IV，因为class3和IV_A共享相同机制）
#'
#' @param dat 任意 generate_dataset_*() 函数的返回值。
#' @export
check_dataset <- function(dat) {
  cat("=== Dataset", dat$dataset, "诊断报告 ===\n")
  cat("性状方差（目标≈1）：", round(apply(dat$Y, 2L, var), 4L), "\n")
  cat("观测遗传力 h²    ：", round(dat$h2_obs, 4L), "\n")

  # 功能位点：有遗传机制但不是工具变量
  n_func <- sum(dat$truth$class %in%
    c("class1", "class2", "class3", "class4", "class5") &
    !dat$truth$is_IV)
  # IV位点：is_IV=TRUE（含标记为class3的IV）
  n_iv <- sum(dat$truth$is_IV)

  cat("功能位点数       ：", n_func, "\n")
  cat("IV辅助位点数     ：", n_iv, "\n")
  cat("Null位点数       ：", sum(dat$truth$class == "null"), "\n")

  # IV强度验证
  if (!is.null(dat$iv_pve_used)) {
    n <- nrow(dat$Y)
    E_F <- dat$iv_pve_used * n / (1 - dat$iv_pve_used)
    cat(sprintf(
      "IV PVE           ：%.4f（期望F统计量 ≈ %.1f，目标>12）\n",
      dat$iv_pve_used, E_F
    ))
  }

  n_show <- min(nrow(dat$truth), 6L + n_iv)
  cat(sprintf("\n真值表（前%d行）：\n", n_show))
  print(dat$truth[seq_len(n_show), ])
  invisible(dat)
}
