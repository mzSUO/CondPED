```
condPED/
│
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── LICENSE
│
├── R/
│   ├── condPED.R                 # 主函数（唯一入口）
│   ├── data_simulation.R          # 数据模拟（复用师姐）
│   ├── qtlnetwork.R              # QTLNetwork 接口

│   ├── projection.R              #  条件投影（OLS实现）
│   ├── gwas.R                    #  GWAS（通用接口）

│   ├── bidirectional.R           # 双向条件 GWAS 流程
│   ├── classification.R          # QTL 分类

│   ├── causal_mr.R               # 因果推断（MR）


│   ├── evaluation.R              # CVR / Gain 等指标
│   ├── visualize.R               # 作图

├── inst/
│   ├── extdata/                  # 示例数据存放用户安装包后读取的辅助数据，小的实例文件
│
├── tests/
│   └── testthat/
	│ ├── test_projection.R 
	│ ├── test_classification.R 
	│ └── test_causal_mr.R
├── man/
│
└── data-raw/
    └── preprocessing_scripts/
    
├── scratch/                # 临时产生的、巨大的中间数据，不占用备份空间
```


# 2. 材料与方法

**Notation**

| **Symbol**                        | **Range**      | **Description**                                                                                                                   |
| --------------------------------- | -------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| $i$                               | $1, \ldots, m$ | Trait index                                                                                                                       |
| $j$                               | $1, \ldots, n$ | Individual index                                                                                                                  |
| $k$                               | $1, \ldots, p$ | Experimental environment index                                                                                                    |
| $l, h$                            | $1, \ldots, s$ | Quantitative Trait SNP (QTS) index                                                                                                |
| $r$                               |                | Number of additive-by-additive epistatic locus pairs                                                                              |
| $c$                               | $1, \ldots, q$ | Covariate index                                                                                                                   |
| $\boldsymbol{\Sigma}_E$           | $m \times m$   | Trait covariance of environmental effects                                                                                         |
| $\boldsymbol{\Sigma}_{AE}$        | $m \times m$   | Trait covariance of additive $\times$ environment effects                                                                         |
| $\boldsymbol{\Sigma}_{AAE}$       | $m \times m$   | Trait covariance of epistasis $\times$ environment effects                                                                        |
| $\boldsymbol{\Sigma}_\varepsilon$ | $m \times m$   | Residual covariance matrix                                                                                                        |
| $\otimes$                         | —              | Kronecker product                                                                                                                 |
| $x_{jl}^a$                        | **基因型变量**      | 个体 $j$ 在位点 $l$ 的**加性基因型编码**。通常采用 0/1/2 编码（即风险等位基因的携带个数）。                                                                          |
| $x_{jlh}^{aa}$                    | **上位性变量**      | 个体 $j$ 在位点 $l$ 与位点 $h$ 之间的**加性 $\times$ 加性上位性编码**。定义为 $x_{jl}^a \times x_{jh}^a$。                                                 |
| $x_{jkc}$                         | **协变量**        | 个体 $j$ 在环境（或群体） $k$ 下第 $c$ 个**协变量**的观测值（如主成分 PC、性别、年龄等）。                                                                          |
| $\theta_{iq}$                     | **通用遗传效应**     | QTS 对性状 $i$ 的**总遗传贡献**。当 $q=l$ 时代表加性效应 $\theta_{il}^a$；当 $q=(l,h)$ 时代表上位性效应 $\theta_{ilh}^{aa}$。                                  |
| $\theta_{iq}^{\text{uncond}}$     | **边际遗传效应**     | **无条件（Marginal）遗传效应**。在不考虑其他性状影响时，通过标准 GWAS 或单表型混合模型估计得到的总效应。                                                                     |
| $\theta_{iq}^{\text{ind}}$        | **独立遗传效应**     | **独立（Independent）遗传效应**。剔除其他性状的介导影响后，QTS 对性状 $i$ 的**直接贡献**。                                                                       |
| $\theta_{iq}^{\text{shared}}$     | **共享遗传效应**     | **共享（Shared）遗传效应**。反映了 QTS 通过性状间的相关性或因果路径对性状 $i$ 产生的间接贡献：$\theta^{\text{shared}} = \theta^{\text{uncond}} - \theta^{\text{ind}}$。 |
|                                   |                |                                                                                                                                   |


## 2.1. Multivariate mixed linear model with additive and epistatic effects

 **Scalar form**

Consider a natural population consists of $n$ individuals and the phenotypic performance of $m$ correlated quantitative traits are measured in $p$ different environments, the genetic variation of the traits are jointly controlled by $s$ quantitative trait QTSs (QTSs) in which $r$ pairs of QTSs are involved in epistatic interactions. To jointly model additive, epistatic, and genotype-by-environment interaction effects, we specify the following mixed linear model:
$$\begin{aligned} y_{ijk} = \mu_i &+ \sum_{c=1}^q b_{ic} x_{jkc} + \sum_{l=1}^s a_{il} x_{jl}^a + \sum_{l=1}^{s-1} \sum_{h=l+1}^{s} aa_{ilh} x_{jlh}^{aa} \\ &+ e_{ik} + \sum_{l=1}^s ae_{ilk} x_{jl}^a + \sum_{l=1}^{s-1} \sum_{h=l+1}^{s} aae_{ilhk} x_{jlh}^{aa} + \varepsilon_{ijk} \end{aligned} \tag{Eq. 1}$$

where:

- $y_{ijk}$ is the phenotypic value of trait $i$ for individual $j$ in environment $k$;
    
- $\mu_i$ is the trait-specific mean;
    
- $b_{ic}$ is the fixed effect of covariate $c$;
    
- $a_{il}$ is the additive effect of locus $l$;
    
- $aa_{ilh}$ is the additive-by-additive epistatic effect between loci $l$ and $h$;
    
- $x_{jl}^a$ and $x_{jlh}^{aa} = x_{jl}^a x_{jh}^a$ denote genotype encodings;
    
- $e_{ik}$, $ae_{ilk}$, and $aae_{ilhk}$ denote environmental, additive-by-environment, and epistasis-by-environment random effects, respectively;$e_{ik} \sim N(0, \sigma_{Ei}^2)，ae_{ilk} \sim N(0, \sigma_{AlE}^2)，aae_{ilhk} \sim N(0, \sigma_{AAlhE}^2)$, respectively
    
- $\varepsilon_{ijk}$ is the residual error，$\varepsilon_{ijk} \sim N(0, \sigma_{\varepsilon i}^2)$. 
    

 All additive and additive-additive epistatic effects are fixed effects, environment effects and interaction of each genetic components with environments are modeled as random effects.


 **Matrix formulation**

Let

$$\mathbf{Y} = \begin{bmatrix} \mathbf{y}_{11} \\ \mathbf{y}_{12} \\ \vdots \\ \mathbf{y}_{np} \end{bmatrix} \in \mathbb{R}^{np \times m},$$

where each row $\mathbf{y}_{jk} \in \mathbb{R}^m$ represents the vector of all traits for individual $j$ in environment $k$.

The model can be written as:

$$\mathbf{Y} = \mathbf{1}\boldsymbol{\mu}^T + \mathbf{X}_c \mathbf{B}_c + \mathbf{X}_A \mathbf{B}_A + \mathbf{X}_{AA} \mathbf{B}_{AA} + \mathbf{U}_E \mathbf{E}_E + \mathbf{U}_{AE} \mathbf{E}_{AE} + \mathbf{U}_{AAE} \mathbf{E}_{AAE} + \mathbf{E}_\varepsilon \tag{Eq. 2}$$

or equivalently:

$$\mathbf{Y} = \mathbf{X}\mathbf{B} + \sum_u \mathbf{U}_u \mathbf{E}_u$$

 **Joint distribution**

To explicitly model cross-trait correlations, the model is formulated in vectorized form:

$$\mathrm{vec}(\mathbf{Y}) \sim N\left( \mathrm{vec}(\mathbf{X}\mathbf{B}), \mathbf{V}_{total} \right)$$

 **Random effect distributions**

All random effects are modeled using matrix normal distributions:

- $\mathbf{E}_E \sim MN(\mathbf{0}, \mathbf{I}_p, \boldsymbol{\Sigma}_E)$
    
- $\mathbf{E}_{AE} \sim MN(\mathbf{0}, \mathbf{I}_{ps}, \boldsymbol{\Sigma}_{AE})$
    
- $\mathbf{E}_{AAE} \sim MN(\mathbf{0}, \mathbf{I}_{pr}, \boldsymbol{\Sigma}_{AAE})$
    
- $\mathbf{E}_\varepsilon \sim MN(\mathbf{0}, \mathbf{I}_{np}, \boldsymbol{\Sigma}_\varepsilon)$
    

where the column covariance matrices $\boldsymbol{\Sigma}_\cdot$ characterize trait–trait covariance structures for each random effect component.

 **Variance–covariance structure**

Under this formulation, the phenotypic variance–covariance matrix is:

$$\mathbf{V}_{total} = \mathbf{U}_E \mathbf{U}_E^\top \otimes \boldsymbol{\Sigma}_E + \mathbf{U}_{AE} \mathbf{U}_{AE}^\top \otimes \boldsymbol{\Sigma}_{AE} + \mathbf{U}_{AAE} \mathbf{U}_{AAE}^\top \otimes \boldsymbol{\Sigma}_{AAE} + \mathbf{I} \otimes \boldsymbol{\Sigma}_\varepsilon $$

**Fixed effects**

The fixed-effect parameter matrix is defined as:

$$\mathbf{B} = (\boldsymbol{\mu}^T, \mathbf{B}_c^T, \mathbf{B}_A^T, \mathbf{B}_{AA}^T)^T$$



## 0.1. 基于条件投影的遗传效应分解框架

Assuming $\mathbf{y}_{jk} \sim N(\boldsymbol{\mu}, \mathbf{V})$, the conditional variance of trait $i$ given the remaining traits $\mathbf{y}_{-i}$ is:

$$
\mathbf{V}_{i \mid -i} = V_{ii} - \mathbf{C}_{i,-i} \mathbf{V}_{-i}^{-1} \mathbf{C}_{-i,i} \tag{Eq. 3}
$$

其中 $\mathbf{C}_{i,-i} = \text{Cov}(y_i, \mathbf{y}_{-i})$ 是性状 $i$ 与其余性状的协方差行向量。
### 0.1.1. 条件表型的构造与性质

While Eq. (3) provides an algebraic expression for the conditional variance $\mathbf{V}_{i \mid -i}$, the conditional variance components (e.g., additive, epistatic) cannot be directly extracted from this matrix form alone (Zhu 1995). Decomposing $\mathbf{V}_{i \mid -i}$ into interpretable genetic components requires fitting a mixed linear model—yet the conditional distribution itself is not directly observable. To overcome this challenge, we adopt an indirect approach: construct an equivalent vector whose variance structure matches $\mathbf{V}_{i \mid -i}$, enabling variance component estimation via standard mixed model machinery.

Estimation of conditional genetic variance components $(\sigma_{u(i | -i)}^{2})$ or prediction of conditional genetic effects $(e_{u(i | -i)})$ cannot be derived directly from the estimate of conditional variance-covariance matrix $(\hat{V}_{(i | -i)})$ . The indirect approaches are suggested for analyzing conditional genetic effects and their variance components.
虽然条件分布本身不可观测，但我们可以构造一个等价的线性投影，其方差结构与真实条件方差完全一致（Eq. 7），从而通过标准混合模型估计条件遗传效应。
Let $\tilde{y}_{ijk} = y_{ijk} - \hat{\mu}_i - \sum_{c=1}^q \hat{b}_{ic}x_{jkc}$ denote the residual phenotype for trait $i$ of individual $j$ in environment $k$, obtained after adjusting for the population mean and covariate effects. We define the conditional phenotype as: $$ \tilde{y}_{ijk}^* = \tilde{y}_{ijk} - \tilde{\mathbf{y}}_{-i,jk}^\top {\gamma}_{i,-i} \tag{Eq. 4} $$ where $\tilde{\mathbf{y}}_{-i,jk} = [\tilde{y}_{1jk},\dots,\tilde{y}_{(i-1)jk},\tilde{y}_{(i+1)jk},\dots,\tilde{y}_{mjk}]^\top$ is the vector of residual phenotypes for all traits excluding $i$,  and the projection coefficient vector is defined as $\gamma_{i,-i} = V_{-i}^{-1} C_{-i,i}$

Here, $V_{-i}$ denotes the phenotypic variance–covariance matrix among traits excluding $i$, and $C_{-i,i}$ is the covariance vector between trait $i$ and the remaining traits.


**Orthogonality property.** From the construction in Eq. (4), the conditional phenotype satisfies the following orthogonality property (see Supplementary Information for proof):
$$
\mathrm{Cov}(\tilde{y}_i^*, \tilde{\mathbf{y}}_{-i}) = \mathbf{0} \tag{Eq. 6}
$$
This result follows from the projection property of the multivariate normal distribution, indicating that the conditional phenotype $\tilde{y}_i^*$ is orthogonal to the remaining traits in the sense of phenotypic covariance.


**Variance equivalence**. The variance of $\tilde{y}_{ijk}^*$ equals the conditional variance in Eq. (3)
$$\text{Var}(\tilde{y}_{ijk}^*) = \mathbf{V}_{ii} - \mathbf{C}_{i,-i} \mathbf{V}_{-i}^{-1} \mathbf{C}_{i,-i} = \mathbf{V}_{i \mid -i} \tag{Eq.7}$$

(derivation in Supplementary Information). Random vector $\tilde{y}_{ijk}^*$ has variance, which is identical to the conditional variance-covariance matrix $\mathbf{V}_{i \mid -i}$ . In practice, unknown parameters in Equation 4 can be replaced by their unbiased estimates.
### 0.1.2. 条件混合模型及遗传效应估计与分解
**Mixed model for conditional phenotypes.** We fit $\tilde{y}_{ijk}^*$ ​ to the same mixed linear model structure as in Eq. (1):
$$\begin{aligned} \tilde{y}_{ijk}^* =  \sum_{l=1}^s a_{il}^* x_{jl}^a + \sum_{l=1}^{s-1} \sum_{h=l+1}^{s} aa_{ilh}^* x_{jlh}^{aa}  \\ + e_{ik}^* + \sum_{l=1}^s ae_{ilk}^* x_{jl}^a + \sum_{l=1}^{s-1} \sum_{h=l+1}^{s}aae_{ilhk}^* x_{jlh}^{aa} + \varepsilon_{ijk}^* \end{aligned} \tag{Eq. 8}$$

 The variance-covariance structure of this model is:

$$\text{Var}(\tilde{y}_{ijk}^*) =  \mathbf{U}_E \mathbf{U}_E^\top \otimes \boldsymbol{\Sigma}_E^* + \mathbf{U}_{AE} \mathbf{U}_{AE}^\top \otimes \boldsymbol{\Sigma}_{AE}^* + \mathbf{U}_{AAE} \mathbf{U}_{AAE}^\top \otimes \boldsymbol{\Sigma}_{AAE}^* + \mathbf{I} \otimes \boldsymbol{\Sigma}_\varepsilon^* =  \mathbf{V}^* \tag{Eq. 9}$$

### 0.1.3. **Distributional equivalence 证明 

Combining Eqs. (7) and (9) yields $\mathbf{V}^* = \mathbf{V}_{i \mid -i}$​. Consequently, variance components estimated from Eq. (8) are unbiased estimators of the conditional variance components:$$\hat{\sigma}_{u^*}^2 = \hat{\sigma}_{u(i \mid -i)}^2 $$

To formally connect the conditional model with the original parameterization,  we define the **conditional random effects** for trait $i$ given the remaining traits as:

$$
e_{(i \mid -i)k}, \quad ae_{(i \mid -i)lk}, \quad aae_{(i \mid -i)lhk}
$$
which represent, respectively, the environmental, additive-by-environment, and epistasis-by-environment effects contributing to the conditional trait variation after removing linear dependence on other traits.
where $u \in \{E, AE, AAE, \varepsilon\}$. Furthermore, the random effects from the conditional model are distributionally equivalent to the true conditional random effects. Formally, let $e_{(i \mid -i)k}$, $ae_{(i \mid -i)lk}$, and $aae_{(i \mid -i)lhk}$ denote the environmental, additive-by-environment, and epistasis-by-environment effects for trait $i$ conditioned on other traits. Then:
$$\begin{aligned} e_{ik}^* &\overset{d}{=} e_{(i \mid -i)k} \sim N(0, \sigma_{E^*}^2), \\ ae_{ilk}^* &\overset{d}{=} ae_{(i \mid -i)lk} \sim N(0, \sigma_{AE^*}^2), \\ aae_{ilhk}^* &\overset{d}{=} aae_{(i \mid -i)lhk} \sim N(0, \sigma_{AAE^*}^2) \end{aligned} \tag{Eq. 10}$$

These conditional random effects can be predicted via best linear unbiased prediction (BLUP; Henderson 1963) or, when variance unbiasedness is required, adjusted unbiased prediction (AUP; Zhu 1995).

**Genetic effect decomposition.** To provide a unified representation of genetic effects, we introduce a generic index $q \in \mathcal{Q}$, where $q = l$ denotes additive effects and $q = (l,h)$ denotes additive-by-additive epistatic effects. Accordingly, all genetic effects are denoted as $\theta_{iq}$.

Under this framework, the conditional phenotype $\tilde{y}_{ijk}^*$ defined in Eq. (4) captures the component of trait $i$ orthogonal to all other trait residuals. Accordingly, any genetic effect estimated from this transformed phenotype represents the component independent of cross-trait linear dependence. We therefore define the independent (orthogonal) genetic effect as: $\theta_{iq}^{\text{ind}} \equiv \theta_{iq}^*$  where $\theta_{iq}^*$ is the effect estimated from the conditional phenotype $\tilde{y}_i^*$. 

By construction, $\theta_{iq}^{\text{ind}}$ represents the trait-specific genetic component that cannot be explained or predicted by other traits. Since the conditional projection represents a linear orthogonal decomposition of the original phenotypic space, the marginal (unconditional) genetic effect can be additively partitioned as: $$\theta_{iq}^{\text{uncond}} = \theta_{iq}^{\text{ind}} + \theta_{iq}^{\text{shared}} \tag{Eq. 11}$$ where $\theta_{iq}^{\text{shared}}$ denotes the shared genetic component that is explained by the linear covariance structure among traits. This decomposition is statistically consistent because the conditional phenotype is constructed to be orthogonal to other traits, ensuring that the two components are statistically uncorrelated and structurally interpretable. Together, they allow us to distinguish between direct, trait-specific genetic effects and indirect, covariance-mediated effects that contribute to pleiotropic genetic architectures.
由于 $\tilde{y}_i^*$ 与 $\tilde{y}_{-i}$ 正交，因此 $\text{Cov}(\theta_{iq}^{\text{ind}}, \theta_{iq}^{\text{shared}}) = 0$，该分解是唯一且无偏的。
## 0.2. Bidirectional conditional GWAS framework


Building upon the conditional projection framework (Section 2.2) and its statistical properties (Section 2.3), we develop a bidirectional conditional testing procedure to identify trait-specific QTL and to decompose their genetic effects into independent and shared components.

1. **A forward conditional analysis**, which removes cross-trait covariance via orthogonal projection to isolate independent genetic effects; and
    
2. **A reverse conditional analysis**, which evaluates whether associations with secondary traits can be fully explained by the target trait, thereby assessing a mediation-like dependency structure.我们通过 reverse conditional test 排除通过其他性状介导的关联

### 0.2.1. Overview of the testing procedure

Let $\mathcal{L}$ denote the set of candidate QTL identified from the multivariate GWAS (Section 2.1). For each trait–locus pair $(i,l)$, the procedure consists of three steps:

1. **Marginal screening (multivariate GWAS)**: detect loci associated with the trait set.
    加入**效应分解**单表型边际效应θiquncond​作为分解的基准，计算独立效应和共享效应逐个性状拟合单表型混合模型
2. **Forward conditional GWAS**: test for trait-specific (independent) effects.
    
3. **Reverse conditional GWAS**: validate whether remaining cross-trait associations are mediated through the target trait.



### 0.2.2. Bidirectional conditional GWAS

**Input:** Phenotype matrix $\mathbf{Y}$, genotype matrix $\mathbf{X}$, covariates $\mathbf{X}_c$, significance thresholds $\alpha_1, \alpha_2, \alpha_3$

**Output:** Classification of QTL–trait pairs into four categories

**Step 1: Marginal multivariate GWAS**

For each locus $l = 1, \ldots, s$, test:

$$H_0: \boldsymbol{\theta}_l = \mathbf{0} \quad \text{vs.} \quad H_1: \boldsymbol{\theta}_l \neq \mathbf{0}$$

using Wilks’ Lambda based on model (2). Define:

$$\mathcal{L} = \{l : P_l^{\mathrm{marg}} < \alpha_1\}, \quad \alpha_1 = QTLNetwork中置换检验的阈值.$$

**Step 2: Forward conditional GWAS **

For each trait $i = 1, \ldots, m$:

1. Construct conditional phenotype:
    
    $$\tilde{y}_{i,jk}^* = \tilde{y}_{i,jk} - \tilde{\mathbf{y}}_{-i,jk}^\top \hat{\mathbf{\gamma}}_{i,-i},$$
    
2. For each $l \in \mathcal{L}$, fit:
    
    $$\tilde{y}_{ijk}^* = \theta_{il}^{\mathrm{ind}} x_{jl} + e_{ik}^* + \varepsilon_{ijk}^*,$$
    

这一步直接在R中编程实现GWAS，注意不是全基因组，是第一步中候选SNP集
    
$\mathcal{S}_1 = \{(l,i): P_{il}^{\mathrm{cond}} < \alpha_2\}, \quad \alpha_2 = \frac{0.05}{|\mathcal{L}| \cdot m}.$
    

**Step 3: Reverse conditional GWAS 

For each each $(l,i) \in \mathcal{S}_1$ and each $t \neq i$, construct the reverse conditional phenotype by conditioning on **all traits except $i$**: $\tilde{y}_{t \mid (-i),jk} = \tilde{y}_{t,jk} - \tilde{\mathbf{y}}_{-i,jk}^\top \boldsymbol{\gamma}_{t,-i}$ 

Construct reverse conditional phenotype:
	    <font color="#ff0000">控制除 i 以外的所有 m−1 个性状后，位点对 t 无效应</font>
> 反向条件检验的核心逻辑是：如果一个 QTS 是性状i特异性的，那么它对任何其他性状t的效应都必须完全通过性状i介导。因此，当我们控制了所有其他性状（除i外）之后，QTS 对t的剩余效应应该完全由i介导，即统计上不显著。
> 
> 反之，如果控制了所有其他性状之后，QTS 对t仍然有显著效应，说明 QTS 对t存在独立于i和所有其他性状的遗传效应，即存在独立多效性。
> 
> 与传统的 " 控制目标性状i" 的方法相比，我们的反向条件设计能够排除其他性状的混杂影响，更准确地识别真正的性状特异性位点。这一设计借鉴了 Byrne 等人 (2020) 在精神疾病遗传学中的研究思路，但扩展到了多环境和上位性效应的情形。

 Fit:
    
    $$ \tilde{y}_{t \mid (-i),jk} = \theta_{t \mid (-i),l} x_{jl} + e_{t \mid (-i),k} + \varepsilon_{t \mid (-i),jk} $$

直接在R中实现GWAS
    
    
 Define Layer 2 significance: $\alpha_3 = \frac{0.05}{|S_1|(m-1)}$.
<font color="#ff0000"> 或者设宽松阈值：For Layer 2 tests, we control the false discovery rate (FDR) at 0.05 across all reverse conditional tests.</font>
A pair $(l,i)$ is declared **trait-specific** if:

$$P_{t \mid (-i),l} \ge \alpha_3, \quad \forall t \neq i$$


### 0.2.3. Interpretation of the bidirectional design

The forward and reverse steps play distinct statistical roles.

- **The forward step** performs a multivariate projection that removes all linear dependence on other traits, yielding an estimate of the independent genetic effect.
    
- The reverse step evaluates whether the association between the locus and secondary traits is statistically eliminated after conditioning on all other traits, which provides a criterion for statistical trait-specificity rather than causal mediation.
    

### 0.2.4. QTL classification

Each QTL–trait pair $(l,i)$ is classified into one of four categories:

| **Category**       | **Criteria**                                                               | **Interpretation**                                                                        |
| ------------------ | -------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| **Not detected**   | $P_{il}^{\mathrm{marg}} \ge \alpha_1$                                      | No significant association with the trait set                                             |
| **Shared**         | $P_{il}^{\mathrm{marg}} < \alpha_1, P_{il}^{\mathrm{cond}} \ge \alpha_2$   | Effect on trait i is fully explained by linear covariance with other traits               |
| **Independent**    | $P_{il}^{\mathrm{cond}} < \alpha_2, \exists t: P_{tl \mid i} < \alpha_3$   | QTS has independent genetic effects on multiple traits                                    |
| **Trait-specific** | $P_{il}^{\mathrm{cond}} < \alpha_2, \forall t: P_{tl \mid i} \ge \alpha_3$ | All genetic effects of the QTS are concentrated on trait iFully concentrated on trait $i$ |
|                    |                                                                            |                                                                             

### 0.3.1. 整合MR的动机

在前述条件投影框架中，我们已经得到：

- SNP 对每个性状的：
    - **边际效应** $\theta^{\text{uncond}}$
    - **条件独立效应** $\theta^{\text{ind}}$

该分解能够区分：

- trait-specific effects
- shared (pleiotropic) effects

但**仍然无法区分：**

- 独立多效性（horizontal pleiotropy）
- 介导结构（vertical / mediated effects）

为此，我们进一步引入**基于个体数据的工具变量回归（MR）**，在不依赖外部 summary data 的情况下，对性状间的**方向性关系进行估计**。

### 0.3.2. 条件独立性作为工具变量选择

令 $n$ 表示样本量，$m$ 表示性状数量，$p$ 表示 SNP 数量。定义：

- $\mathbf{G} \in \mathbb{R}^{n \times p}$：标准化基因型矩阵（每列均值为 0，方差为 1）
    
- $\mathbf{y}_i \in \mathbb{R}^n$：所有个体中性状 $i$ 的向量
    
- $\mathbf{Y} = [\mathbf{y}_1, \ldots, \mathbf{y}_m] \in \mathbb{R}^{n \times m}$：表型矩阵
    

对于假设存在 $A \to B$ 因果关系的性状对 $(A, B)$：

- 性状 $A$ 是 **暴露 (exposure)**
    
- 性状 $B$ 是 **结果 (outcome)**
    
- 我们的目标是估计因果效应 $\gamma_{A \to B}$
    

#### 0.3.2.1. 边缘效应与条件效应（简述）
In the causal inference analysis, we focus on additive genetic effects  ($\theta^a$), which serve as valid instrumental variables under the  Mendelian randomization framework.
根据第 2.2 节，SNP $l$ 对性状 $i$ 的边缘效应为：

$$\hat{\theta}_{il}^{\text{uncond}} = \frac{\mathbf{g}_l^T \mathbf{y}_i}{\mathbf{g}_l^T \mathbf{g}_l} = \frac{1}{n} \mathbf{g}_l^T \mathbf{y}_i$$

条件表型（公式 4）产生独立效应：

$$\hat{\theta}_{il}^{\text{ind}} = \frac{1}{n} \mathbf{g}_l^T \mathbf{y}_{i|-i}$$

其中 $\mathbf{y}_{i|-i} = \mathbf{y}_i - \mathbf{Y}_{-i} \hat{\boldsymbol{\gamma}}_{i,-i}$，且 $\hat{\boldsymbol{\gamma}}_{i,-i} = (\mathbf{Y}_{-i}^T \mathbf{Y}_{-i})^{-1} \mathbf{Y}_{-i}^T \mathbf{y}_i$。

---

### 0.3.3. 通过条件独立性选择工具变量

传统的孟德尔随机化要求工具变量满足：

1. **相关性 (Relevance)**：SNP 与暴露相关 ($G_l \to A$)
    
2. **独立性 (Independence)**：无混杂因素影响 ($G_l \perp\!\!\!\perp U$)
    
3. **排他性约束 (Exclusion restriction)**：对结果无直接影响（$G_l$ 仅通过 $A$ 影响 $B$）
    

标准 MR 方法通过 $P$ 值测试相关性，但无法直接评估排他性约束。我们使用条件效应来解决这个问题。

#### 0.3.3.1. 选择标准

如果 SNP $l$ 满足以下条件，则被选为 $A \to B$ 的有效工具变量：

$$\mathcal{Z}_A = \left\{ l : \begin{array}{l} P_{Al}^{\text{marg}} < \alpha_{\text{rel}} \quad \text{(相关性)} \\[0.5ex] P_{B|A,l}^{\text{cond}} > \alpha_{\text{excl}} \quad \text{(排他性)} \\[0.5ex] |r_{ll'}| < r_{\text{LD}} \ \forall l' \in \mathcal{Z}_A, l' < l \quad \text{(LD 聚合)} \end{array} \right\}$$

其中：

- $P_{Al}^{\text{marg}}$ 是检验 $H_0: \theta_{Al}^{\text{uncond}} = 0$ 的 $P$ 值。
    
- $P_{B|A,l}^{\text{cond}}$ 是检验 $H_0: \theta_{B|A,l}^{\text{ind}} = 0$ 的 $P$ 值（来自条件 GWAS）。
    
- 典型阈值：$\alpha_{\text{rel}} = 5 \times 10^{-8}$，$\alpha_{\text{excl}} = 0.05$，$r_{\text{LD}} = 0.1$。
    

**核心创新**：排他性约束通过 $\theta_{B|A,l}^{\text{ind}}$ 进行 **直接** 检验，它量化了在移除对 $A$ 的依赖后 SNP 对 $B$ 的效应。如果 $\theta_{B|A,l}^{\text{ind}} \approx 0$，则该 SNP 在构造上满足排他性约束。

令 $k = |\mathcal{Z}_A|$ 表示有效工具变量的数量，$\mathbf{G}_A \in \mathbb{R}^{n \times k}$ 为相应的基因型子矩阵。

### 0.3.4. 2.5.4.2SLS估计与双向检验

使用多工具变量时，我们采用两阶段最小二乘法 (2SLS)：

#### 0.3.4.1. 第一阶段：预测暴露

$$\mathbf{y}_A = \mathbf{G}_A \boldsymbol{\alpha}_A + \boldsymbol{\varepsilon}_A$$

OLS 估计值为：

$$\hat{\boldsymbol{\alpha}}_A = (\mathbf{G}_A^T \mathbf{G}_A)^{-1} \mathbf{G}_A^T \mathbf{y}_A$$

拟合值为：

$$\hat{\mathbf{y}}_A = \mathbf{G}_A \hat{\boldsymbol{\alpha}}_A = \mathbf{P}_A \mathbf{y}_A$$

其中 $\mathbf{P}_A = \mathbf{G}_A (\mathbf{G}_A^T \mathbf{G}_A)^{-1} \mathbf{G}_A^T$ 是向 $\mathbf{G}_A$ 列空间投射的投影矩阵。

#### 0.3.4.2. 第二阶段：估计因果效应

$$\mathbf{y}_B = \gamma_{A \to B} \hat{\mathbf{y}}_A + \boldsymbol{\varepsilon}_B$$

2SLS 估计量为：

$$\hat{\gamma}_{A \to B}^{\text{2SLS}} = \frac{\mathbf{y}_B^T \mathbf{P}_A \mathbf{y}_A}{\mathbf{y}_A^T \mathbf{P}_A \mathbf{y}_A} = (\mathbf{y}_A^T \mathbf{P}_A \mathbf{y}_A)^{-1} \mathbf{y}_A^T \mathbf{P}_A \mathbf{y}_B$$

---

### 0.3.5. 使用汇总统计量重新表述


#### 0.3.5.1. SNP 级别效应向量

定义：

$$\hat{\boldsymbol{\theta}}_A = [\hat{\theta}_{Al_1}^{\text{marg}}, \ldots, \hat{\theta}_{Al_k}^{\text{marg}}]^T \in \mathbb{R}^k$$

$$\hat{\boldsymbol{\theta}}_B = [\hat{\theta}_{Bl_1}^{\text{marg}}, \ldots, \hat{\theta}_{Bl_k}^{\text{marg}}]^T \in \mathbb{R}^k$$

其中 $l_1, \ldots, l_k \in \mathcal{Z}_A$ 是选定的工具变量。

#### 0.3.5.2. LD 相关矩阵

样本内 LD 结构由下式捕捉：

$$\mathbf{R}_A = \frac{1}{n} \mathbf{G}_A^T \mathbf{G}_A \in \mathbb{R}^{k \times k}$$

其中 $R_{A,ij} = \text{cor}(\mathbf{g}_{l_i}, \mathbf{g}_{l_j})$ 是 SNP $l_i$ 和 $l_j$ 之间的 Pearson 相关系数。

#### 0.3.5.3. 等价于加权汇总统计量形式

将 $\hat{\theta}_{il}^{\text{uncond}}$ 和 $\mathbf{P}_A$ 的定义代入 2SLS 公式，并利用 $\mathbf{y}_i^T \mathbf{g}_l = n \hat{\theta}_{il}^{\text{uncond}}$，我们得到：

$$\hat{\gamma}_{A \to B}^{\text{2SLS}} = \frac{\hat{\boldsymbol{\theta}}_A^T \mathbf{R}_A^{-1} \hat{\boldsymbol{\theta}}_B}{\hat{\boldsymbol{\theta}}_A^T \mathbf{R}_A^{-1} \hat{\boldsymbol{\theta}}_A}$$

这等价于权重为 $\mathbf{W} = \mathbf{R}_A^{-1}$ 的广义最小二乘法 (GLS)，反映了工具变量之间由 LD 引起的决策结构。

---

### 0.3.6. 方差估计

在标准 2SLS 渐近线下，$\hat{\gamma}_{A \to B}^{\text{2SLS}}$ 的方差为：

$$\text{Var}(\hat{\gamma}_{A \to B}^{\text{2SLS}}) = \frac{\hat{\sigma}_B^2}{\hat{\boldsymbol{\theta}}_A^T \mathbf{R}_A^{-1} \hat{\boldsymbol{\theta}}_A}$$



$$\hat{\sigma}_B^2 = \frac{1}{n-k-1} (\mathbf{y}_B - \hat{\gamma}_{A \to B}^{\text{2SLS}} \hat{\mathbf{y}}_A)^T (\mathbf{y}_B - \hat{\gamma}_{A \to B}^{\text{2SLS}} \hat{\mathbf{y}}_A)$$

#### 0.3.6.1. 为什么 LD 很重要

LD 矩阵 $\mathbf{R}_A$ 捕捉了工具变量之间的相关性。忽略 LD（即使用 $\mathbf{W} = \mathbf{I}$）会导致标准误估计偏倚和推断失效。

---

### 0.3.7. 异质性检验

为了评估工具变量是否提供一致的估计值，我们计算 Cochran's $Q$ 统计量。定义单个 Wald 比率：

$$\hat{\theta}_i = \frac{\hat{\theta}_{Bl_i}^{\text{marg}}}{\hat{\theta}_{Al_i}^{\text{marg}}}, \quad i = 1, \ldots, k$$

则：

$$Q = \sum_{i=1}^k \frac{(\hat{\theta}_i - \hat{\gamma}_{A \to B}^{\text{2SLS}})^2}{\text{SE}(\hat{\theta}_i)^2} \sim \chi^2_{k-1} \quad (\text{在 } H_0 \text{ 下})$$

较大的 $Q$ 值 ($P < 0.05$) 可能暗示存在残余多效性、非线性关系或群体分层。此时建议进行敏感性分析。

---

### 0.3.8. 双向估计

对于性状对 $(A, B)$，我们可以测试两个方向：

- $A \to B$：使用工具变量 $\mathcal{Z}_A$ 估计 $\hat{\gamma}_{A \to B}^{\text{2SLS}}$
    
- $B \to A$：使用工具变量 $\mathcal{Z}_B$ 估计 $\hat{\gamma}_{B \to A}^{\text{2SLS}}$
    

显著性的不对称性 ($T_{A \to B} \gg T_{B \to A}$) 为方向性依赖提供了证据。
## 2.5. 性状间的因果推断

### 2.5.1. 整合MR的动机

在前述条件投影框架中，我们已经得到：

- SNP 对每个性状的：
    - **边际效应** $\theta^{\text{uncond}}$
    - **条件独立效应** $\theta^{\text{ind}}$

该分解能够区分：

- trait-specific effects
- shared (pleiotropic) effects

但**仍然无法区分：**

- 独立多效性（horizontal pleiotropy）
- 介导结构（vertical / mediated effects）

为此，我们进一步引入**基于个体数据的工具变量回归（MR）**，在不依赖外部 summary data 的情况下，对性状间的**方向性关系进行估计**。

### 2.5.2. 通过条件独立性选择工具变量

在QTLNetwork的混合线性模型中，QTS $l$ 对性状 $i$ 的边际效应 $\hat{\theta}_{il}^{\text{uncond}}$ 通过REML/MINQUE估计获得：

传统的孟德尔随机化要求工具变量满足：

1. **相关性 (Relevance)**：SNP 与暴露相关 ($G_l \to A$)
    
2. **独立性 (Independence)**：无混杂因素影响 ($G_l \perp\!\!\!\perp U$)
    
3. **排他性约束 (Exclusion restriction)**：对结果无直接影响（$G_l$ 仅通过 $A$ 影响 $B$）
    

标准 MR 方法通过 $P$ 值测试相关性，但无法直接评估排他性约束。我们使用条件效应来解决这个问题。


如果 SNP $l$ 满足以下条件，则被选为 $A \to B$ 的有效工具变量：

$$\mathcal{Z}_A = \left\{ l : \begin{array}{l} P_{Al}^{\text{marg}} < \alpha_{\text{rel}} \quad \text{(相关性)} \\[0.5ex] P_{B|A,l}^{\text{cond}} > \alpha_{\text{excl}} \quad \text{(排他性)} \\[0.5ex] |r_{ll'}| < r_{\text{LD}} \ \forall l' \in \mathcal{Z}_A, l' < l \quad \text{(LD 聚合)} \end{array} \right\}$$

其中：

- $P_{Al}^{\text{marg}}$ 是检验 $H_0: \theta_{Al}^{\text{uncond}} = 0$ 的 $P$ 值。
    
- $P_{B|A,l}^{\text{cond}}$ 是检验 $H_0: \theta_{B|A,l}^{\text{ind}} = 0$ 的 $P$ 值（来自条件 GWAS）。
    
- 典型阈值：$\alpha_{\text{rel}} = 5 \times 10^{-8}$，$\alpha_{\text{excl}} = 0.05$，$r_{\text{LD}} = 0.1$。
    

**核心创新**：排他性约束通过 $\theta_{B|A,l}^{\text{ind}}$ 进行 **直接** 检验，它量化了在移除对 $A$ 的依赖后 SNP 对 $B$ 的效应。如果 $\theta_{B|A,l}^{\text{ind}} \approx 0$，则该 SNP 在构造上满足排他性约束。

令 $k = |\mathcal{Z}_A|$ 表示有效工具变量的数量，$\mathbf{G}_A \in \mathbb{R}^{n \times k}$ 为相应的基因型子矩阵。

### 2.5.3. LD-aware因果效应估计

设：

$\hat{\boldsymbol{\theta}}_A,\ \hat{\boldsymbol{\theta}}_B$ 分别为暴露性状A与结果性状B在工具变量集合上的效应向量。

我们通过如下GLS问题估计因果效应：

$$\hat{\gamma}_{A \to B} = \arg\min_{\gamma} (\hat{\boldsymbol{\theta}}_B - \gamma \hat{\boldsymbol{\theta}}_A)^T \mathbf{R}^{-1} (\hat{\boldsymbol{\theta}}_B - \gamma \hat{\boldsymbol{\theta}}_A)$$

其解析解为：

$$\hat{\gamma}_{A \to B} = \frac{\hat{\boldsymbol{\theta}}_A^T \mathbf{R}^{-1} \hat{\boldsymbol{\theta}}_B}{\hat{\boldsymbol{\theta}}_A^T \mathbf{R}^{-1} \hat{\boldsymbol{\theta}}_A}$$
### 2.5.4. 双向估计

对于性状对 $(A, B)$，我们可以测试两个方向：

- $A \to B$：使用工具变量 $\mathcal{Z}_A$ 估计 $\hat{\gamma}_{A \to B}^{\text{2SLS}}$
    
- $B \to A$：使用工具变量 $\mathcal{Z}_B$ 估计 $\hat{\gamma}_{B \to A}^{\text{2SLS}}$
    

显著性的不对称性 ($T_{A \to B} \gg T_{B \to A}$) 为方向性依赖提供了证据。



## 2.7. Model evaluation and interpretation metrics

To facilitate biological interpretation of the proposed decomposition framework, we introduce a set of summary statistics that quantify the relative contributions of independent and shared genetic effects, as well as trait-level genetic architecture.

**Contribution ratio**

The relative contribution of independent genetic effects is quantified as:

$$\rho_{il}^{\mathrm{ind}} = \frac{(\theta_{iq}^{\mathrm{ind}})^2}{(\theta_{iq}^{\mathrm{ind}})^2 + (\theta_{iq}^{\mathrm{shared}})^2}$$

which measures the proportion of variance attributable to trait-specific effects.

**Conditional variance ratio (CVR)**

At the trait level, we define:

$$\mathrm{CVR}_{i \mid -i} = \frac{\sigma_{i \mid -i}^2}{\sigma_i^2}$$

which quantifies the proportion of genetic variation that is independent of other traits. A larger CVR indicates a stronger trait-specific genetic basis.

**Detection gain**

To evaluate the practical benefit of conditional analysis, we define:

$$\mathrm{Gain}_i = \frac{|\{l: P_{il}^{\mathrm{cond}} < \alpha_2, \, P_{il}^{\mathrm{marg}} \ge \alpha_1\}|}{|\{l: P_{il}^{\mathrm{marg}} < \alpha_1\}|}$$

which measures the proportion of loci detectable only after removing cross-trait confounding.

## 2.8. Simulation study 

为了系统评估所提出的 **CondPED** (_Conditional Projection-based Effect Decomposition_) 方法的统计性能，我们设计了一系列模拟实验，从以下四个方面进行验证：

- **检测能力 (Power)**：评估方法识别性状特异性遗传效应的能力。
- **I 类错误控制 (Type I error)**：检验在零假设下的错误率控制。
- **遗传效应分解准确性 (Decomposition accuracy)**：评估独立效应与共享效应的估计精度。
- **多效性结构识别能力 (Pleiotropy classification)**：区分不同遗传机制（独立效应、共享效应、介导效应、连锁伪多效性）。

### 模拟场景

| **场景** | **描述**                  | **生物学目的**                                                       | **统计学目的**                            | **数学设定**                                                                                                                                                                                                                                                                                                                                       | **预期结果**                                                                                              |
| ------ | ----------------------- | --------------------------------------------------------------- | ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- |
| **S0** | **Null scenario**       | N/A                                                             | 验证 Type I error control              | 所有 $\boldsymbol{\beta}_l = \mathbf{0}$                                                                                                                                                                                                                                                                                                         | Type I error $\le 0.05$<br>                                                                           |
| **S1** | **Trait-specific QTL**  | 性状特异性基因<br>                                                     | 验证条件投影不削弱真实信号  <br>检验独立效应检测能力        | $\boldsymbol{\beta}_l = (0, \ldots, \beta_{li}, \ldots, 0)$<br>只有第 $i$ 个元素非零                                                                                                                                                                                                                                                                   | Power 高<br><br>正确分类为 "trait_specific"<br><br><br>$\beta_l^{\text{ind}} \approx \beta_l^{\text{true}}$ |
| **S2** | **Fully shared QTL**    | 泛效性基因<br><br>                                                   | 验证条件投影能消除共享成分                        | $\boldsymbol{\beta}_l = (\beta, \beta, \ldots, \beta)$<br><br>潜变量显式构造共享结构                                                                                                                                                                                                                                                                      | 条件 GWAS 后信号消失<br>分类为 "shared"<br>$\beta_l^{\text{ind}} \approx \mathbf{0}$                            |
| **S3** | **Mediated pleiotropy** | 介导效应<br><br>  <br><br>(如：SNP $\rightarrow$ 肥胖 $\rightarrow$ 血压) | 区分直接/间接效应<br>验证 Forward + Reverse 框架 | **完全介导 (S3a)：**<br>$y_1 = x_l \beta_1 + e_1$<br>$y_2 = \gamma { x_l \beta_1} + e_2$<br>$\Rightarrow \beta_{l2}^{\text{marg}} = \gamma\beta_1, \beta_{l2}^{\text{ind}} = 0$<br>**部分介导 (S3b)：**<br>$y_2 = \gamma y_1 + x_l \beta_2 + e_2$<br>$\Rightarrow \beta_{l2}^{\text{marg}} = \gamma\beta_1 + \beta_2, \beta_{l2}^{\text{ind}} = \beta_2$ | **S3a:** Trait 1 为 "trait_specific"<br>Trait 2 为 "shared"<br>**S3b:** 两者都有独立成分<br>分类为 "independent"   |
| **S4** | **LD-linked loci**      | 连锁不平衡<br><br>(假多效性)                                             | 区分真实多效性 vs LD 伪多效性                   | $\text{SNP}_1 \to \text{Trait}_1 (\beta_{11} \neq 0, \beta_{12} = 0)$<br><br>$\text{SNP}_2 \to \text{Trait}_2 (\beta_{21} = 0, \beta_{22} \neq 0)$<br>$\text{corr}(\text{SNP}_1, \text{SNP}_2) = r_{LD}$ (如 0.6)- 一个 LD block（5–10 SNP）<br>- causal SNP 在 block 中                                                                              | 两个 SNP 都被分类为 "trait_specific"<br><br>不会被误判为真实多效性                                                      |



| **场景**     | **样本量 (n)**       | **特征数 (m)** | **相关性 (ρ)**        | **效应大小**        | **重复次数** | **实验目的**    |
| ---------- | ----------------- | ----------- | ------------------ | --------------- | -------- | ----------- |
| **主分析**    | 1000              | 4           | 0.6                | 0.3             | 1000     | 展示各场景分类准确性  |
| **样本量敏感性** | {500, 1000, 2000} | 4           | 0.6                | 0.3             | 500      | 展示 Power 曲线 |
| **相关性敏感性** | 1000              | 4           | {0, 0.3, 0.6, 0.9} | 0.3             | 500      | 展示方法稳健性     |
| **高维挑战**   | 1000              | {2, 4, 8}   | 0.6                | 0.3             | 500      | 展示可扩展性      |
| **效应大小**   | 1000              | 4           | 0.6                | {0.1, 0.3, 0.5} | 500      | 展示检测阈值      |



### 检测指标


| 指标                           |                                                                                        |
| ---------------------------- | -------------------------------------------------------------------------------------- |
| power-检测能力                   | $\text{Power} = \frac{\text{正确检测到的 QTL 数}}{\text{真实 QTL 数}}$                           |
| FDR                          |                                                                                        |
| I 类错误率 (Type I error)        | $\text{Type I error} = \Pr(P < \alpha)$                                                |
| Classification accuracy（最重要） | $\text{Accuracy} = \frac{\text{正确分类的 QTL 数}}{\text{总 QTL 数}}$                          |
| trait-specific识别             |                                                                                        |
| shared vs independent区分      |                                                                                        |
| 遗传效应分解误差                     | $\text{MSE} = \|\hat{\theta}^{\mathrm{ind}} - \theta^{\mathrm{ind}}_{\text{true}}\|^2$ |
| 检测增益 (Detection gain)        | $\mathrm{Gain}_i$                                                                      |

