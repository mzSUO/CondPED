2.1. Multivariate mixed linear model with additive and epistatic effects

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
    

 All additive and additive-additive epistatic effects are fixed effects, environment effects and interaction of each genetic components with environments are modeled as random effects. 模型的矩阵形式及方差成分估计细节详见补充材料S1。
## 0.1. 基于条件投影的遗传效应分解框架

一个边际SNP–性状关联并不一定代表该位点对目标性状存在真实的独立遗传作用，该关联还可能来源于性状间的协方差传播，或反映位点同时作用于多个性状的多效性结构。此外，多效性本身又可能对应不同的生物学机制，包括位点对多个性状的独立直接作用（horizontal pleiotropy），或通过性状间因果路径传递产生的间接作用（vertical pleiotropy）。

为区分上述不同情形，本研究将条件投影分析与孟德尔随机化（Mendelian randomization, MR）整合于统一的混合模型框架中。首先通过条件投影将边际遗传效应正交分解为独立遗传效应与协方差介导效应；随后通过跨性状独立效应比较识别单性状特异性位点与多效性位点；最后利用 MR 推断多效性位点的因果作用机制。整个框架实现了从统计混杂到直接效应，再到因果机制的层级化解耦。

### 0.1.1. 理论框架

令 $\mathbf{y}_{jk} \in \mathbb{R}^m$ 表示个体 $j$ 在环境 $k$ 下所有性状的观测向量。在混合模型框架下，性状间的边际表型协方差矩阵为：

$$\mathbf{V} = \boldsymbol{\Sigma}_E + \boldsymbol{\Sigma}_{AE} + \boldsymbol{\Sigma}_{AAE} + \boldsymbol{\Sigma}_\varepsilon$$

该矩阵聚合了所有随机效应在性状层面的协方差成分。根据多元正态分布的条件分布理论，性状 $i$ 在给定其余性状 $\mathbf{y}_{-i}$ 条件下的条件方差为：

$$\mathbf{V}_{i \mid -i} = V_{ii} - \mathbf{C}_{i,-i} \mathbf{V}_{-i}^{-1} \mathbf{C}_{-i,i}$$

其中 $\mathbf{C}_{i,-i}$ 为性状 $i$ 与其余性状的协方差向量。该条件方差对应于$\mathbf{V}_{-i}$ 在表型协方差矩阵中的Schur补，在几何上表示移除了与其他性状线性依赖后的残余变异。

### 0.1.2. 条件表型的构造与性质

While Eq. (3) provides an algebraic expression for the conditional variance $\mathbf{V}_{i \mid -i}$, the conditional variance components (e.g., additive, epistatic) cannot be directly extracted from this matrix form alone (Zhu 1995). Decomposing $\mathbf{V}_{i \mid -i}$ into interpretable genetic components requires fitting a mixed linear model—yet the conditional distribution itself is not directly observable. To overcome this challenge, we adopt an indirect approach: construct an equivalent vector whose variance structure matches $\mathbf{V}_{i \mid -i}$, enabling variance component estimation via standard mixed model machinery.

> Estimation of conditional genetic variance components $(\sigma_{u(i | -i)}^{2})$ or prediction of conditional genetic effects $(e_{u(i | -i)})$ cannot be derived directly from the estimate of conditional variance-covariance matrix $(\hat{V}_{(i | -i)})$ . The indirect approaches are suggested for analyzing conditional genetic effects and their variance components.

Let $\tilde{y}_{ijk} = y_{ijk} - \hat{\mu}_i - \sum_{c=1}^q \hat{b}_{ic}x_{jkc}$ denote the residual phenotype for trait $i$ of individual $j$ in environment $k$, obtained after adjusting for the population mean and covariate effects. We define the conditional phenotype as: $$ \tilde{y}_{ijk}^* = \tilde{y}_{ijk} - \tilde{\mathbf{y}}_{-i,jk}^\top {\gamma}_{i,-i} \tag{Eq. 4} $$ where $\tilde{\mathbf{y}}_{-i,jk} = [\tilde{y}_{1jk},\dots,\tilde{y}_{(i-1)jk},\tilde{y}_{(i+1)jk},\dots,\tilde{y}_{mjk}]^\top$ is the vector of residual phenotypes for all traits excluding $i$,  and the projection coefficient vector is defined as $\gamma_{i,-i} = V_{-i}^{-1} C_{-i,i}$. 该构造确保 $\tilde{y}_i^*$ 与其余性状在表型协方差意义上正交，即 $\mathrm{Cov}(\tilde{y}_i^*, \tilde{\mathbf{y}}_{-i}) = \mathbf{0}$，且满足方差等价性 $\mathrm{Var}(\tilde{y}_i^*) = \mathbf{V}_{i \mid -i}$(证明见补充材料)

对条件表型拟合与边际模型相同的混合线性模型结构，记条件模型中的遗传效应为 $\theta_{iq}^*$（其中 $q=l$ 表示加性效应，$q=(l,h)$ 表示上位性效应）。根据补充材料中证明的方差等价性，我们构造的条件表型的方差严格等于给定所有其他性状后的理论表型条件方差$\mathbf{V}_{i \mid -i}$ ，从条件模型估计的方差成分 $\hat{\sigma}_{u^*}^2$ 是真实条件方差成分 $\sigma_{u(i \mid -i)}^2$ 的无偏估计（$u \in \{E, AE, AAE, \varepsilon\}$）。进一步地，条件模型中的随机效应与真实条件随机效应在分布上等价（补充材料）。

**Mixed model for conditional phenotypes.** We fit $\tilde{y}_{ijk}^*$ ​ to the same mixed linear model structure as in Eq. (1):
$$\begin{aligned} \tilde{y}_{ijk}^* =  \sum_{l=1}^s a_{il}^* x_{jl}^a + \sum_{l=1}^{s-1} \sum_{h=l+1}^{s} aa_{ilh}^* x_{jlh}^{aa}  \\ + e_{ik}^* + \sum_{l=1}^s ae_{ilk}^* x_{jl}^a + \sum_{l=1}^{s-1} \sum_{h=l+1}^{s}aae_{ilhk}^* x_{jlh}^{aa} + \varepsilon_{ijk}^* \end{aligned} \tag{Eq. 8}$$

 The variance-covariance structure of this model is:

$$\text{Var}(\tilde{y}_{ijk}^*) =  \mathbf{U}_E \mathbf{U}_E^\top \otimes \boldsymbol{\Sigma}_E^* + \mathbf{U}_{AE} \mathbf{U}_{AE}^\top \otimes \boldsymbol{\Sigma}_{AE}^* + \mathbf{U}_{AAE} \mathbf{U}_{AAE}^\top \otimes \boldsymbol{\Sigma}_{AAE}^* + \mathbf{I} \otimes \boldsymbol{\Sigma}_\varepsilon^* =  \mathbf{V}^* \tag{Eq. 9}$$
由此得到条件遗传效应。

### 0.1.3. 遗传效应的分解
由于条件表型在表型协方差意义上与所有其他性状正交，从其估计的遗传效应代表了**不依赖于其他性状线性预测的成分**。我们将从条件表型估计的遗传效应定义为**条件遗传效应** $\theta_{iq}^{\mathrm{cond}} \equiv \theta_{iq}^*$。基于条件投影的线性性质，边际遗传效应可唯一分解为两个统计上不相关的成分：。


$$\theta_{iq}^{\mathrm{marg}} = \theta_{iq}^{\mathrm{cond}} + \theta_{iq}^{\mathrm{cov}}$$

其中 $\theta_{iq}^{\mathrm{cov}}$ 表示由性状间线性协方差结构解释的共享成分。由于 $\tilde{y}_i^*$ 与 $\tilde{\mathbf{y}}_{-i}$ 正交，两个成分在统计上不相关（$\mathrm{Cov}(\theta_{iq}^{\mathrm{cond}}, \theta_{iq}^{\mathrm{cov}}) = 0$）。This decomposition represents a statistical orthogonal reparameterization of SNP effects induced by phenotype-space projection, rather than a decomposition of underlying biological causal mechanisms.
## 0.2. 分层遗传效应解耦

基于上述条件投影框架，我们建立了一个统一的分层检测流程，旨在系统性地将多性状GWAS中的边际关联分解为可解释的生物学模式。

### 0.2.1. 检测流程设计

令 $\mathcal{L}$ 表示经置换检验（1,000次）确定的全基因组显著位点集合。对于每个候选对 $(l,i)$（$l \in \mathcal{L}$, $i=1,\ldots,m$），依次执行四个步骤：

**第一层：边际多性状GWAS关联检测**。利用QTLNetwork中的多性状混合模型检验 $H_0: \boldsymbol{\theta}_l = \mathbf{0}$，识别显著位点集合 $\mathcal{L} = \{l: P_l^{\mathrm{marg}} < \alpha_1\}$，其中 $\alpha_1$ 为经验全基因组阈值。对每个性状 $i$ 进一步拟合单性状模型，估计边际效应 $\hat{\theta}_{il}^{\mathrm{marg}}$ 作为后续分解基准。
借鉴 Byrne 等人 (2020) 在精神疾病遗传学中的研究策略，我们仅对原始边际 GWAS 中已达到全基因组显著性的位点进行后续条件分析。仅比较边际模型与条件模型的效应值差异不足以可靠识别具有统计显著性的遗传效应，因为对相关性状具有强效应的位点在条件化后可能出现统计上的伪效应变化。通过明确评估条件效应的统计显著性，我们避免了对条件效应位移的误读。

>

**第二层：前向条件分析独立遗传效应识别**。对每个性状 $i$，按2.2.2节所述构造条件表型 $\tilde{y}_{i,jk}^*$，并对 候选位点$l \in \mathcal{L}$ 拟合条件混合模型，估计独立效应 $\hat{\theta}_{il}^{\mathrm{cond}}$。采用Bonferroni校正的显著性阈值 $\alpha_2 = 0.05/(|\mathcal{L}| \times m)$。我们将候选位点初步划分为两类：
- 协方差解释关联：位点对性状i的边际效应显著，但独立效应不显著 ($P_{il}^{\mathrm{marg}} < \alpha_1, P_{il}^{\mathrm{ind}} \ge \alpha_2$)，表明该关联可由性状间协方差完全解释，no detectable covariance-independent genetic effect； 
- 条件显著关联：位点对性状i的独立效应显著 ($P_{il}^{\mathrm{ind}} < \alpha_2$)，表明该关联不能被性状间协方差完全解释，位点存在不能被其他性状线性解释的遗传效应。covariance-explained component

**第三层：单性状与多性状效应区分**

对于具有条件显著关联的位点，我们通过比较其在所有性状上的条件遗传效应的统计显著性，进一步区分单性状特异性模式与多效性模式：
- **性状特异性条件效应**：位点仅在一个性状上具有显著的条件遗传效应 ($P_{il}^{\mathrm{ind}} < \alpha_2$ 且 $\forall t \neq i, P_{tl}^{\mathrm{ind}} \ge \alpha_2$)，表明该位点仅在目标性状上检测到与其他性状正交的遗传成分，在所有其他性状上均未检测到显著的条件遗传效应；
- **多效性效应模式**：位点在两个或更多性状上具有显著的条件遗传效应 ($\exists t \neq i, P_{il}^{\mathrm{ind}} < \alpha_2$ 且 $P_{tl}^{\mathrm{ind}} < \alpha_2$)，表明该位点在多个性状上均检测到与其他性状正交的遗传成分。

**第四层：多效性的因果机制分解**

多效性位点可能对应两种不同的生物学机制：**水平多效性**，即位点独立直接影响多个性状，不存在性状间的因果依赖；以及**垂直多效性**，即位点仅直接影响上游性状，其对下游性状的效应通过性状间的因果通路介导传递。为区分这两种机制，我们整合孟德尔随机化方法进行因果推断。

**工具变量筛选**

传统孟德尔随机化的核心挑战在于难以验证排他性约束，即工具变量仅通过暴露性状影响结果性状，而对结果性状无独立于暴露的直接效应。我们利用条件 GWAS 的结果为这一假设提供经验证据：对于因果方向 $A \rightarrow B$，若位点对 B 的独立遗传效应 $\theta_{B \mid A,l}^{\mathrm{ind}}$ 统计上不显著，则与孟德尔随机化的排他性约束一致，为 "该位点对 B 没有独立于 A 的直接效应" 这一假设提供了经验支持。
（我们的结果没有违反排他性约束，但是不能说证明了排他性约束，不要过度声明）

据此，我们定义有效工具变量集合为：
$$\mathcal{Z}_A = \left\{ l : \begin{array}{l} P_{Al}^{\mathrm{marg}} < 5 \times 10^{-8} \quad \text{(相关性)} \\ P_{B \mid A,l}^{\mathrm{ind}} > 0.05 \quad \text{(排他性)} \\ |r_{ll'}| < 0.1 \ \forall l' \in \mathcal{Z}_A, l' < l \quad \text{(LD修剪)} \end{array} \right\}$$
需要特别强调的是，MR 分析使用的是全基因组独立的显著位点，而非当前正在分类的单个位点，避免了循环推断问题。同时，条件效应不显著仅表明在当前样本量和统计模型下未检测到显著的直接效应，不能绝对证明水平多效性不存在。

**考虑 LD 的广义最小二乘因果估计**

为校正工具变量间残留连锁不平衡导致的估计偏倚，我们采用广义最小二乘法 (GLS) 估计因果效应。令 $\hat{\boldsymbol{\theta}}_A$ 和 $\hat{\boldsymbol{\theta}}_B$ 分别为工具变量集合对暴露性状 A 和结果性状 B 的边际效应向量，$\mathbf{R}$ 为工具变量间的 LD 相关矩阵。因果效应的 GLS 估计为：
$$\hat{\gamma}_{A \rightarrow B} = \frac{\hat{\boldsymbol{\theta}}_A^\top \mathbf{R}^{-1} \hat{\boldsymbol{\theta}}_B}{\hat{\boldsymbol{\theta}}_A^\top \mathbf{R}^{-1} \hat{\boldsymbol{\theta}}_A} \tag{7}$$
采用 Wald 检验评估因果效应的统计显著性。对于每一对性状，都进行了双向多因素分析（bidirectional multivariate analysis）。如果某个基因位点的因果效应在某个方向上显着（例如，A → B，P < 0.05），但在相反方向上不显着（B → A，P ≥ 0.05），则该基因位点被归类为具有“垂直多效性”（vertical pleiotropy）。如果两个方向的效应都显着，则该基因位点被标记为“双向多效性”或“存在混杂因素”的类型，并被排除在“水平多效性”（horizontal pleiotropy）或“垂直多效性”的分类之外；因为这种模式可能反映了相互影响的关系或未测量的混杂因素。如果两个方向的效应都不显着，则该基因位点被归类为具有“水平多效性”。

### 0.2.2. 介导机制的一致性验证

对于 MR 支持的因果通路 $A \rightarrow B$，我们进一步进行反向条件分析以验证介导机制的一致性：

$$\tilde{y}_{t \mid i,jk} = \tilde{y}_{t,jk} - \gamma_{ti} \tilde{y}_{i,jk} \tag{8}$$

在经典的完全介导模型，控制上游性状 $i$ 后，位点对下游性状 $t$ 的效应应该显著减弱或消失。这一分析为因果推断提供了额外的支持证据。

对于MR支持的因果通路，进行了逆向条件分析以验证介导一致性。我们没有采用任意比例阈值，而是通过单侧Wilcoxon符号秩检验（检验配对差异（原始效应减去调整后效应）是否随机大于零，评估调整上游性状后，遗传效应对下游性状的效应是否系统性减弱。P值<0.05被视为显著效应衰减的证据，因此通过了一致性验证。
### 0.2.3. 统一分类体系

综合上述四层分析结果，我们将每个性状-位点对归入以下五个类别：

**表 1. SNP-性状关联的分层分类体系**

| **类别**                                             | **判定标准统计条件**                                                                                                                | **生物学解释**                                                                                                                            |                                                   |
| -------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------- |
| **未检测到**                                           | $P_{il}^{\mathrm{marg}} \ge \alpha_1$ 对所有性状 $i$                                                                             | 边际 GWAS 未达显著性阈值。该位点与任何性状均无显著遗传关联，或效应低于检测阈值                                                                                           | 位点对所有性状无作用，或作用过弱                                  |
| **协方差介导关联**covariance-induced /explained component | $P_{il}^{\mathrm{marg}} < \alpha_1$ 且 $P_{il}^{\mathrm{ind}} \ge \alpha_2$                                                  | 边际显著但独立效应不显著。位点仅对其他性状有真实遗传效应，目标性状上的显著关联完全是性状间遗传协方差被动牵连的统计假象                                                                          | 该位点与目标性状的边际关联可完全由其他性状的线性协方差解释，未检测到与其他性状正交的遗传成分    |
| **性状特异性独立效应**                                      | $P_{il}^{\mathrm{ind}} < \alpha_2$ 且 $\forall t \neq i, P_{tl}^{\mathrm{ind}} \ge \alpha_2$                                 | 仅在一个性状上具有显著独立效应。位点的独立遗传效应完全集中在单一性状上，对其他性状无独立作用                                                                                       | 该位点仅在一个性状上检测到显著的条件遗传效应，在所有其他性状上均未检测到显著的正交遗传成分     |
| **水平多效性**                                          | $\exists t \neq i, P_{il}^{\mathrm{ind}} < \alpha_2$ 且 $P_{tl}^{\mathrm{ind}} < \alpha_2$，且 MR 不支持任何性状间的因果关系                | 多个性状具有显著独立效应，但无因果证据。位点独立且直接地影响多个性状，不存在性状间的上下游因果关系，反映基因的多功能性                                                                          | 该位点在多个性状上检测到显著的条件遗传效应，且无统计证据支持性状间存在因果介导关系         |
| **垂直多效性**                                          | $\exists t \neq i, P_{il}^{\mathrm{ind}} < \alpha_2$ 且 $P_{tl}^{\mathrm{ind}} < \alpha_2$，且 MR 显著支持 $i \rightarrow t$ 的因果关系 | 多个性状具有显著独立效应，且存在因果证据。位点对下游性状的效应<font color="#ff0000">存在</font>通过上游性状的因果通路介导传递，揭示了遗传变异作用的上下游传递机制，不能说<font color="#ff0000">完全中介</font> | 该位点在多个性状上检测到显著的条件遗传效应，且统计证据与 "上游性状因果介导下游性状" 的模型一致 |
| 双向多效性                                              | 两个方向都显著                                                                                                                     |                                                                                                                                      |                                                   |
## 0.3. 统计检验与计算实现


{全基因组关联分析}。边际多性状 GWAS 采用 QTLNetwork 2.0 软件实现，该软件基于混合线性模型框架联合检验多个性状的遗传效应。为控制实验整体 I 类错误率，临界 F 值通过 Henderson III 方法结合置换检验（1,000 次重复）在 0.05 显著性水平下确定。对所有显著信号进行逐步回归筛选，排除虚假关联并识别条件独立的 QTS 集合 $\mathcal{L}$，用于后续条件分析。

{条件分析与因果推断}。所有后续分析在 R 4.4.2 环境下实现。表型协方差矩阵 $\mathbf{V}$ 及投影系数 $\boldsymbol{\gamma}_{i,-i}$ 通过残差表型的样本协方差估计。方差成分采用 REML 算法估计。显著 QTS 的遗传效应通过 
**解释性指标如下表**。可视化使用包实现。

模型评估指标体系

| 评估维度      | 指标                          | 数学定义                                                                                                                                  | 解释                                            |
| :-------- | :-------------------------- | :------------------------------------------------------------------------------------------------------------------------------------ | :-------------------------------------------- |
| **统计有效性** | I类错误率 $\alpha_{\text{emp}}$ | $\frac{\sum_{l \in \mathcal{H}_0} \mathbb{1}(P_l < \alpha)}{\|\mathcal{H}_0\|}$                                                       | 零假设位点的假阳性率，目标 ≤ 0.05                          |
|           | 统计效能 Power                  | $\frac{\sum_{l \in \mathcal{H}_1} \mathbb{1}(P_l < \alpha)}{\|\mathcal{H}_1\|}$                                                       | 真阳性位点的检出率，目标 ≥ 0.80                           |
|           | 检测增益 Gain$_i$               | $\frac{\|{l: P_{il}^{\mathrm{cond}} < \alpha_2, P_{il}^{\mathrm{marg}} \geq \alpha_1}\|}{\|{l: P_{il}^{\mathrm{marg}} < \alpha_1}\|}$ | 条件分析中显著但边际分析中不显著的位点比例，反映了条件分析在控制协方差混杂后的检测能力提升 |
| **分类准确性** | 整体准确率 Acc                   | $\frac{1}{N} \sum_{(l,i)} \mathbb{1}(\hat{C}_{li} = C_{li}^{\text{true}})$                                                            | 五分类判定全部正确的比例                                  |
|           | 类别F1分数 $F1_c$               | $\frac{2 \cdot \text{Prec}_c \cdot \text{Rec}_c}{\text{Prec}_c + \text{Rec}_c}$                                                       | 每个类别的精确率与召回率调和平均                              |
| **估计精度**  | 效应估计RMSE                    | $\sqrt{\frac{1}{N} \sum_{(l,i)} (\hat{\theta}_{il}^{\mathrm{cond}} - \theta_{il}^{\text{true}})^2}$                                   | 条件遗传效应估计的均方根误差                                |

**符号说明**：$\mathcal{H}_0$ 为零假设位点集，$\mathcal{H}_1$ 为模拟设定的效应位点集，$C_{li}$ 为QTL-性状对的模拟设定的分类标签分类标签（未检测到 / 协方差解释关联 / 性状特异性条件效应 / 水平多效性模式 / 推定垂直多效性），$\mathbb{1}(\cdot)$ 为指示函数。