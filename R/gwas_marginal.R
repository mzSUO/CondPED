# ==============================================================================
# gwas_marginal.R — Layer 1：多性状边际关联扫描
# ==============================================================================
#
# 【统计模型】（论文 Eq. 1，仅加性固定效应，OLS）
#
#   y_{ij} = mu_i + sum_c b_{ic} * x_{jc} + a_{il} * x_{jl} + eps_{ij}
#
#   eps_{ij} ~ N(0, sigma²_{eps,i})
#
# 【设计理念：为什么 Layer 1 用 OLS 而不是混合模型】
#   Layer 1 是筛选步骤，不是最终推断步骤。使用 OLS 的原因：
#   (1) 单环境数据：随机效应（e_{ik}, ae_{ilk}）不可识别，混合模型会模型错设。
#   (2) 多环境数据：环境作为固定协变量（哑变量编码）控制主效应，
#       不需要估计方差成分（方差成分属于 Layer 2 推断范畴）。
#   (3) 速度与透明度：OLS 可处理 p=100,000 个 SNP；混合模型 REML 慢10-100倍，
#       在筛选阶段不提供额外价值。
#
# 【有意排除的效应项】
#   上位性和 G×E 随机效应（aa_{ilh}, ae_{ilk}, aae_{ilhk}，对应论文 Eq.1）
#   有意排除在 Layer 1 之外。这些高阶效应仅在 Layer 2 的条件模型中估计，
#   且只针对 |L|<<p 的小候选集，避免多重检验负担。
#
# 【多重检验校正】
#   Bonferroni 校正，覆盖所有（SNP × 性状）对：alpha / (p × m)
#   这控制了"所有加性效应均为零"这一联合零假设的家族错误率（FWER）
#
# 【输出】
#   - 每个（SNP，性状）对的效应估计、标准误、p值
#   - 位点级别摘要：n_marg（有显著信号的性状数）
#   - 分类：class1（n_marg=1）vs. layer2候选（n_marg≥2）
# ==============================================================================


# ==============================================================================
# 第一部分：内部 OLS 计算引擎（向量化，一次处理所有 SNP）
# ==============================================================================

#' 单性状向量化单SNP-OLS扫描（内部函数）
#'
#' 对 X 的每一列同时拟合 y ~ 1 + Z + x，使用闭合OLS公式，
#' 避免逐SNP调用 lm() 的开销。速度约为循环版本的50-100倍。
#'
#' 【算法步骤】
#'   1. QR分解偏出协变量（含截距），得到残差化的 y_r 和 X_r
#'   2. 向量化 OLS：beta_l = (x_r_l' y_r) / (x_r_l' x_r_l)
#'   3. 利用 RSS_l = ||y_r||² - beta_l² × (x_r_l' x_r_l) 计算残差方差
#'   4. 计算 SE 和 t 统计量，导出双侧 p 值
#'
#' 【偏出协变量的含义】
#'   与每个SNP模型都包含协变量等价，但只做一次QR分解（计算量从 p 次降为1次）。
#'   QR残差化等价于 Frisch-Waugh-Lovell 定理的应用。
#'
#' @param y  数值向量，长度n。单个性状的表型值。
#' @param X  n×p 数值矩阵。基因型矩阵（加性编码）。
#' @param Z  n×q 数值矩阵，协变量矩阵（不含截距）。
#'   NULL=仅包含截距。多环境时传入哑变量矩阵。
#' @return 数据框，p行，包含列：beta（效应量）、se（标准误）、
#'   t_stat（t统计量）、p_value（双侧p值）。
#' @keywords internal
.ols_scan_one_trait <- function(y, X, Z = NULL) {
    n <- nrow(X)
    p <- ncol(X)

    # ---- 步骤1：通过QR分解偏出协变量（含截距）----------------------------------
    # W 包含截距列（全1）和额外协变量（如环境哑变量）
    # 只做一次 QR 分解，等价于对所有 p 个 SNP 都控制了协变量
    if (is.null(Z)) {
        W <- matrix(1, n, 1) # 仅截距
    } else {
        W <- cbind(1, Z) # 截距 + 协变量（环境哑变量等）
    }

    qrW <- qr(W)
    y_r <- qr.resid(qrW, y) # 残差化表型：去掉协变量解释的部分
    X_r <- qr.resid(qrW, X) # 残差化基因型：n×p 矩阵（一次性处理所有SNP）

    # 每个 SNP 模型的残差自由度 = n - (协变量数 + 截距) - 1（SNP本身）
    df_res <- n - ncol(W) - 1L

    # ---- 步骤2：向量化 OLS 估计：beta = (X_r'y_r) / diag(X_r'X_r) -----------
    # 利用矩阵运算同时估计所有 p 个 SNP 的效应量
    XtX_diag <- colSums(X_r^2) # 长度p，每个SNP的 x_r'x_r（标量）
    Xty <- crossprod(X_r, y_r) # p×1，每个SNP的 x_r'y_r

    beta <- as.numeric(Xty) / XtX_diag # 向量化除法，长度p

    # ---- 步骤3：计算残差方差和标准误 ------------------------------------------
    # RSS_l = ||y_r||² - beta_l² × (x_r_l'x_r_l)
    # 避免构造 n×p 的拟合值矩阵（内存高效）
    y_r_ss <- sum(y_r^2) # 标量，与l无关
    rss <- y_r_ss - beta^2 * XtX_diag # 长度p，每个SNP的残差平方和
    rss <- pmax(rss, 0) # 数值保护（防止浮点负值）

    sigma2 <- rss / df_res # 残差方差，长度p
    se <- sqrt(sigma2 / XtX_diag) # beta的标准误，长度p

    # ---- 步骤4：t统计量和双侧p值 -----------------------------------------------
    t_stat <- beta / se
    p_value <- 2 * pt(abs(t_stat), df = df_res, lower.tail = FALSE)

    data.frame(
        beta = beta,
        se = se,
        t_stat = t_stat,
        p_value = p_value,
        stringsAsFactors = FALSE
    )
}


# ==============================================================================
# 第二部分：Layer 1 主函数
# ==============================================================================

#' Layer 1 — 多性状边际关联扫描
#'
#' 对每个SNP、每个性状执行 OLS 扫描（论文 Eq.1，仅加性固定效应），
#' 采用 Bonferroni 校正（alpha/(p×m)）控制所有（SNP，性状）对的家族错误率。
#' 根据每个位点在多少个性状上显著，将候选位点分为两类：
#'   - class1（n_marg=1）：仅一个性状显著，直接标注为性状特异性关联
#'   - layer2（n_marg≥2）：多性状显著，作为多效性候选位点进入 Layer 2
#'
#' 【多环境支持】
#'   传入 env 参数时，自动将环境转为哑变量编码并加入协变量矩阵 Z。
#'   环境作为固定效应，吸收不同环境的整体表型偏移，保持 OLS 框架。
#'   这与 Layer 2 的混合模型（随机截距+随机斜率）不同：
#'   Layer 1 只需控制主效应（固定），无需估计方差成分（随机）。
#'
#' 【MAF过滤】
#'   扫描前过滤低频 SNP（默认 MAF<0.01），避免极端低频变异导致数值不稳定。
#'   MAF 计算方式兼容 RIL（-1/1编码）和 F2（-1/0/1编码）。
#'
#' @param X       n×p 数值矩阵。基因型矩阵，需有列名（SNP ID）。
#'   支持 RIL（-1/1）、F2（-1/0/1）等加性编码方式。
#' @param Y       n×m 数值矩阵。表型矩阵，需有列名（性状名）。
#'   个体顺序必须与 X 一致。
#' @param Z       n×q 数值矩阵，额外固定协变量（不含截距），或 NULL。
#'   若同时提供 env，环境哑变量会自动追加到 Z 后面。
#' @param env     因子或字符向量，长度n，指示每个观测所属的环境。
#'   NULL=单环境数据（默认）。
#'   提供时：自动哑变量编码（丢弃第一个水平避免共线性），
#'   环境主效应作为固定协变量纳入，等价于论文 Eq.1 中的 e_{ik} 固定部分。
#' @param alpha   家族显著性水平，默认0.05。
#'   Bonferroni 阈值 = alpha / (p × m)。
#' @param min_maf MAF 过滤阈值，默认0.01。
#'   低于此阈值的 SNP 在扫描前被移除。
#' @param verbose 逻辑值。是否打印进度信息，默认TRUE。
#'
#' @return 命名列表，包含：
#' \describe{
#'   \item{class1}{字符向量。仅在一个性状上显著（n_marg=1）的SNP ID。
#'     直接标注为 Class 1（性状特异性关联），不进入 Layer 2。}
#'   \item{layer2}{字符向量。在两个及以上性状上显著（n_marg≥2）的SNP ID。
#'     作为多效性候选位点，进入 Layer 2 条件分析。}
#'   \item{snp_summary}{数据框，每个显著SNP一行。
#'     列：snp_id、n_marg、以及每个性状的 beta_i、se_i、p_i、sig_i。
#'     按 n_marg 降序排列（多效性候选在前）。}
#'   \item{scan_result}{长格式数据框，每个（SNP，性状）对一行，
#'     包含所有SNP（含不显著）的效应估计和原始p值。
#'     用于 Layer 3 MR 的工具变量筛选（需要全量扫描结果）。}
#'   \item{threshold}{数值。实际使用的 Bonferroni 校正阈值。}
#'   \item{n_tested}{整数。通过 MAF 过滤后实际扫描的 SNP 数。}
#'   \item{alpha}{数值。传入的显著性水平 alpha。}
#' }
#'
#' @export
layer1_marginal_scan <- function(X, Y,
                                 Z = NULL,
                                 env = NULL,
                                 alpha = 0.05,
                                 min_maf = 0.01,
                                 verbose = TRUE) {
    # ---- 0. 输入验证 -----------------------------------------------------------
    stopifnot(
        is.matrix(X), is.numeric(X),
        is.matrix(Y), is.numeric(Y),
        nrow(X) == nrow(Y), # 个体数必须一致
        !is.null(colnames(X)), # SNP ID 必须存在
        !is.null(colnames(Y)) # 性状名必须存在
    )

    n <- nrow(X)
    p <- ncol(X)
    m <- ncol(Y)
    snp_ids <- colnames(X)
    trait_ids <- colnames(Y)

    if (verbose) {
        message(sprintf(
            "[Layer 1] n=%d | p=%d SNPs | m=%d 性状 | alpha=%.3f",
            n, p, m, alpha
        ))
    }

    # ---- 1. 构建协变量矩阵 Z（处理环境效应）-----------------------------------
    # 截距在 .ols_scan_one_trait 内部通过 QR 步骤处理，这里只追加额外协变量

    if (!is.null(env)) {
        env <- factor(env)

        if (nlevels(env) < 2L) {
            # 只有一个环境水平，哑变量无法构造，退化为单环境模式
            warning("[Layer 1] env 只有一个水平，忽略环境协变量。")
            env <- NULL
        } else {
            # 哑变量编码：丢弃第一个水平（避免与截距共线）
            # K个环境 → K-1个哑变量列
            env_dummies <- model.matrix(~env)[, -1, drop = FALSE]
            colnames(env_dummies) <- paste0("env_", levels(env)[-1])

            # 追加到已有协变量矩阵后面（或直接作为Z）
            Z <- if (is.null(Z)) env_dummies else cbind(Z, env_dummies)

            if (verbose) {
                message(sprintf(
                    "[Layer 1] 检测到%d个环境，已添加%d个哑变量协变量。",
                    nlevels(env), ncol(env_dummies)
                ))
            }
        }
    }

    # ---- 2. MAF 过滤 -----------------------------------------------------------
    # 计算每个 SNP 的次等位基因频率（兼容 RIL/-1/1 和 F2/-1/0/1 编码）
    col_means <- colMeans(X)
    # 将编码映射到[0,1]频率空间（max-min标准化）
    allele_freq <- (col_means - min(X)) / (max(X) - min(X))
    maf_vec <- pmin(allele_freq, 1 - allele_freq) # 取较小的等位基因频率

    keep <- maf_vec >= min_maf # 逻辑向量：TRUE=保留
    n_remove <- sum(!keep)
    if (n_remove > 0 && verbose) {
        message(sprintf("[Layer 1] 移除 %d 个 MAF<%.3f 的 SNP。", n_remove, min_maf))
    }

    X_scan <- X[, keep, drop = FALSE] # 过滤后的基因型矩阵
    snp_keep <- snp_ids[keep] # 保留的 SNP ID
    p_scan <- sum(keep) # 实际扫描的 SNP 数

    # ---- 3. Bonferroni 显著性阈值 ---------------------------------------------
    # 校正因子 = p×m（所有SNP×性状对的数量）
    # 这控制了所有对联合零假设（全部无效应）的 FWER
    threshold <- alpha / (p_scan * m)

    if (verbose) {
        message(sprintf(
            "[Layer 1] Bonferroni 阈值：%.3e = %.3f / (%d × %d)",
            threshold, alpha, p_scan, m
        ))
    }

    # ---- 4. 逐性状 OLS 扫描 ---------------------------------------------------
    # 对每个性状调用一次向量化 OLS，同时处理所有 p 个 SNP
    # 结果存入 p×m 矩阵（beta、se、p值）
    mat_beta <- matrix(NA_real_, p_scan, m, dimnames = list(snp_keep, trait_ids))
    mat_se <- matrix(NA_real_, p_scan, m, dimnames = list(snp_keep, trait_ids))
    mat_pval <- matrix(NA_real_, p_scan, m, dimnames = list(snp_keep, trait_ids))

    for (j in seq_len(m)) {
        if (verbose) message(sprintf("  扫描性状 %d/%d：%s ...", j, m, trait_ids[j]))
        res_j <- .ols_scan_one_trait(y = Y[, j], X = X_scan, Z = Z)
        mat_beta[, j] <- res_j$beta
        mat_se[, j] <- res_j$se
        mat_pval[, j] <- res_j$p_value
    }

    # ---- 5. 显著性矩阵 --------------------------------------------------------
    # mat_sig[l, i] = TRUE 表示 SNP l 在性状 i 上显著（p < Bonferroni阈值）
    mat_sig <- mat_pval < threshold # 逻辑矩阵，p×m

    # 每个 SNP 显著关联的性状数（n_marg）
    n_marg <- rowSums(mat_sig) # 长度p_scan

    # ---- 6. 位点分类 ----------------------------------------------------------
    # 只有至少一个性状显著的 SNP 进入候选集
    sig_snps <- snp_keep[n_marg >= 1L]

    # Class 1 候选：仅在一个性状上显著
    # 可能是真性状特异性位点，也可能是完全中介的Class 3（需Layer 3 MR区分）
    class1 <- snp_keep[n_marg == 1L]

    # Layer 2 候选：在两个及以上性状上显著
    # 进入条件分析（Layer 2）区分水平多效/垂直多效/双向因果
    layer2 <- snp_keep[n_marg >= 2L]

    if (verbose) {
        message(sprintf(
            "[Layer 1] 显著位点：%d | Class 1：%d | Layer 2候选：%d",
            length(sig_snps), length(class1), length(layer2)
        ))
    }

    # ---- 7. 组装宽格式摘要表（snp_summary）-------------------------------------
    # 每个显著 SNP 一行，包含所有性状的效应、SE、p值和显著性标记
    if (length(sig_snps) > 0L) {
        base_df <- data.frame(
            snp_id = sig_snps,
            n_marg = n_marg[n_marg >= 1L],
            stringsAsFactors = FALSE
        )

        # 按性状追加列（宽格式）：beta_Trait1, se_Trait1, p_Trait1, sig_Trait1, ...
        for (j in seq_len(m)) {
            tid <- trait_ids[j]
            base_df[[paste0("beta_", tid)]] <- mat_beta[sig_snps, j]
            base_df[[paste0("se_", tid)]] <- mat_se[sig_snps, j]
            base_df[[paste0("p_", tid)]] <- mat_pval[sig_snps, j]
            base_df[[paste0("sig_", tid)]] <- mat_sig[sig_snps, j]
        }

        # 按 n_marg 降序排列（多性状多效性位点排在前面，便于快速查看）
        snp_summary <- base_df[order(base_df$n_marg, decreasing = TRUE), ]
        rownames(snp_summary) <- NULL
    } else {
        # 无显著位点时返回空数据框并给出警告
        snp_summary <- data.frame(snp_id = character(0), n_marg = integer(0))
        warning("[Layer 1] 未检出显著位点，考虑放宽 alpha 或增大样本量。")
    }

    # ---- 8. 组装长格式全量扫描结果（scan_result）-------------------------------
    # 包含所有 SNP（含不显著）的所有性状结果
    # 用途：
    #   - Layer 3 MR 工具变量筛选（相关性准则需要全量p值）
    #   - 阈值敏感性分析
    #   - QQ图绘制（Null位点p值）
    scan_result <- data.frame(
        snp_id = rep(snp_keep, times = m), # 每个SNP重复m次
        trait = rep(trait_ids, each = p_scan), # 每个性状重复p次
        beta = as.vector(mat_beta), # 按列展开：先SNP1的m个性状，再SNP2...
        se = as.vector(mat_se),
        p_value = as.vector(mat_pval),
        sig = as.vector(mat_sig),
        stringsAsFactors = FALSE
    )

    # ---- 9. 返回结果列表 -------------------------------------------------------
    list(
        class1      = class1, # Class 1候选（n_marg=1）
        layer2      = layer2, # Layer 2候选（n_marg≥2）
        snp_summary = snp_summary, # 宽格式摘要表
        scan_result = scan_result, # 长格式全量结果
        threshold   = threshold, # Bonferroni阈值
        n_tested    = p_scan, # 实际扫描SNP数
        alpha       = alpha # 输入的显著性水平
    )
}


# ==============================================================================
# 第三部分：验证辅助函数
# ==============================================================================

#' 对照模拟真值验证 Layer 1 结果
#'
#' 根据模拟数据集的 truth 表，计算 Layer 1 的 Power 和 FDR。
#'
#' 【Power 定义】
#'   Power = 被检出的真实功能位点数 / 真实功能位点总数
#'   真实功能位点 = truth$class %in% true_classes（如 class1, class3, class4 等）
#'
#' 【FDR 定义】
#'   FDR = 被错误检出的位点数 / 总检出位点数
#'   错误检出 = 检出位点不在真实功能位点集合中（注意：IV位点有真实效应，
#'   不应计为误检——simulate_data.R 已将其标注为 class3）
#'
#' @param l1           layer1_marginal_scan() 的输出结果。
#' @param truth        模拟数据的真值表（generate_dataset_*()$truth），
#'   需有 SNP 和 class 列。
#' @param true_classes 字符向量。哪些 class 标签被视为真实功能位点。
#'   默认包含全部五类：class1 到 class5。
#' @return 不可见返回列表，包含 power、fdr、per_locus（逐位点检出情况）。
#' @export
validate_layer1 <- function(l1,
                            truth,
                            true_classes = c(
                                "class1", "class2",
                                "class3", "class4", "class5"
                            )) {
    true_qtl <- truth$SNP[truth$class %in% true_classes] # 所有真实功能位点
    detected <- c(l1$class1, l1$layer2) # 所有Layer 1检出位点

    tp <- sum(true_qtl %in% detected) # 真阳性
    fn <- sum(!true_qtl %in% detected) # 假阴性（漏检）
    fp <- sum(!detected %in% true_qtl) # 假阳性（误检）
    pwr <- tp / length(true_qtl) # Power
    fdr <- if (length(detected) > 0) fp / length(detected) else NA_real_ # FDR

    cat(sprintf(
        "\n=== Layer 1 验证报告 ===\n\n  真实功能位点数 : %d\n  检出位点数     : %d\n  真阳性         : %d / %d  (Power = %.1f%%)\n  假阳性         : %d       (FDR   = %.1f%%)\n  漏检           : %s\n\n",
        length(true_qtl), length(detected),
        tp, length(true_qtl), pwr * 100,
        fp, fdr * 100,
        paste(setdiff(true_qtl, detected), collapse = ", ")
    ))

    # 逐位点检出情况表（用于详细诊断）
    per_locus <- data.frame(
        SNP = true_qtl,
        true_class = truth$class[match(true_qtl, truth$SNP)],
        detected = true_qtl %in% detected,
        stringsAsFactors = FALSE
    )

    invisible(list(power = pwr, fdr = fdr, per_locus = per_locus))
}
