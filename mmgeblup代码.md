```
#' calculateKaXa function
#'
#' @param geno_data genotype data
#' @param pheno_data phenotype data
#' @param A additive relationship matrix
#' @param site_qtl a vector, names of QTL sites
#' @param m_qtl the number of QTLs
#'
#' @return list(Ka = Ka,Xa = Xa)
#' @export
#' @import dplyr
#'
#' @examples \dontrun{KaXa = calculateKaXa(mmgeno_data, mmpheno_data, A, site_qtl, m_qtl)}
calculateKaXa <- function(geno_data, pheno_data, A, site_qtl, m_qtl)
{
  # extract qtl information from what we got from GWAS
  if(m_qtl>0)
  {
    Ka = A.mat(as.matrix(geno_data[,!colnames(geno_data) %in% site_qtl]))
    Xa = geno_data[as.character(pheno_data$GID), site_qtl]
    Xa = as.matrix(Xa)
    colnames(Xa) = paste0(site_qtl,"_A")
  } else {
    Ka = A
    Xa = c()
  }

  colnames(Ka) = rownames(Ka) = rownames(geno_data)

  return(list(Ka = Ka,Xa = Xa))
}
```


```
#' calculateKdXd function
#'
#' @param geno_data genotype data
#' @param pheno_data phenotype data
#' @param D dominance relationship matrix
#' @param site_qtl_dom a vector, names of dom QTL sites
#' @param m_qtl_dom the number of dom QTLs
#'
#' @return list(Kd = Kd,Xd = Xd)
#' @export
#' @import dplyr
#'
#' @examples \dontrun{KdXd = calculateKdXd(mmgeno_data, mmpheno_data, D, site_qtl, m_qtl)}
calculateKdXd <- function(geno_data, pheno_data, D, site_qtl_dom, m_qtl_dom)
{
  # extract qtl information from what we got from GWAS
  if(m_qtl_dom>0)
  {
    Kd = D.mat(as.matrix(geno_data[,!colnames(geno_data) %in% site_qtl_dom]))
    Xd = 1-abs(geno_data[as.character(pheno_data$GID), site_qtl_dom]) # transform -1,0,1 to 0,1,0 code
    Xd = as.matrix(Xd)
    colnames(Xd) = paste0(site_qtl_dom,"_D")
  } else {
    Kd = D
    Xd = c()
  }

  colnames(Kd) = rownames(Kd) = rownames(geno_data)

  return(list(Kd = Kd,Xd = Xd))
}
```

```
#' gblup function
#'
#' @param data a data frame
#' @param A additive genetic relationship matrix
#' @param D dominance genetic relationship matrix
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = gblup(model = "AD", data = dt, A = mmdata$A, D = mmdata$D)}
gblup <- function(model, data, A, D)
{
  ## Evaluate input
  if (missing(A)) {
    stop("Error: no input for additive kinship matrix")
  }

  if (model == "A") {
    mod = mmer(reformulate("1", "trait"),
               random = ~vsr(GID, Gu = A),
               rcov = ~units,
               data = data,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    BV = subset.data.frame(data, select = "GID")
    BV$mu = mod$Beta$Estimate[1]  # mu
    BV$A = mod$U$`u:GID`$trait[BV$GID]  # A
    BV$pre = BV$mu + BV$A
  } else if (model == "AD") {
    if (missing(D)) {
      stop("Error: no input for dominance kinship matrix")
    }

    datafake = data %>% dplyr::mutate(GID1=GID)

    mod = mmer(reformulate("1", "trait"),
               random = ~vsr(GID, Gu = A) + vsr(GID1, Gu = D),
               rcov = ~units,
               data = datafake,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    BV = subset.data.frame(data, select = "GID")
    BV$mu = mod$Beta$Estimate[1]  # mu
    BV$A = mod$U$`u:GID`$trait[BV$GID]  # A
    BV$D = mod$U$`u:GID1`$trait[BV$GID]  # D
    BV$pre = BV$mu + BV$A + BV$D
  } else {
    stop("Error: Model type not recognized. Please provide either 'A' or 'AD'.")
  }

  return(list(mod, BV))
}
```


```
#' geno.generate function
#' A simple function to generate homozygotic genotype data (coded with 1 or -1)
#' with respect to the number of individual, the number of markers, MAF range and chromosome.
#' Markers are independent.
#'
#' @param indNum The number of individual
#' @param snpNum The number of markers
#' @param maf.min The lower bound of MAF, so that MAF is sampled from U(maf.min, maf.max)
#' @param maf.max The upper bound of MAF, so that MAF is sampled from U(maf.min, maf.max)
#' @param chr.snpNum The number of snps on each chromosome
#'
#' @return geno_data a data frame, where the first three columns are $CHR, $SNP, and $BP
#' @export
#' @importFrom stats rbinom runif
#'
#' @examples geno.generate(100, 200, 0.05, 0.5, c(50,50,40,60))
geno.generate <- function(indNum, snpNum, maf.min, maf.max, chr.snpNum)
{

  if(sum(chr.snpNum)!=snpNum){
    stop("Error: geno.generate() function: the summation of chr.snpNum is not equal to snpNum.")
  }

  Ga = matrix(NA, indNum, snpNum)
  freq = runif(snpNum, maf.min, maf.max)
  for(j in 1:snpNum)
  {
    Ga[,j] = rbinom(indNum, 2, freq[j])-1              # additive(-1,0,1)
    # Ga[,j] = ifelse(rbinom(indNum, 1, freq[j])==0,-1,1)  # additive(-1,1)
  }
  rownames(Ga) = paste0("gid", 1:indNum)
  colnames(Ga) = paste0("snp", 1:snpNum)

  chrName = paste0("chr", 1:length(chr.snpNum))


  geno_data = data.frame(CHR = as.vector(unlist(sapply(chrName, function(x) rep(x, chr.snpNum[which(chrName == x)])))),
                         SNP = colnames(Ga),
                         BP = as.vector(unlist(sapply(chr.snpNum, function(x) seq(1:x)))),
                         t(Ga))

  return(geno_data)
}
```


```
#' mmdata function
#' Format data for mmGEBLUP analysis
#'
#' @param geno_data input genotype data
#' @param pheno_data input phenotype data
#' @param qtl_data input additice qtl data
#' @param qtl_dom_data input dominance qtl data
#'
#' @return a list with QCed genotype and phenotype data, and other data required for GS analysis
#'         $ mmgeno_data an n*m matrix
#'         $ mmpheno_data an n_obs*3 data frame
#'         $ Ka additive relationship matrix for mmGBLUP
#'         $ Xa additive fixed effect coefficient for mmGBLUP
#'         $ Kd dominance relationship matrix for mmGBLUP
#'         $ Xd dominance fixed effect coefficient for mmGBLUP
#'         $ A additive relationship matrix for GBLUP
#'         $ D additive relatinoship matrix for GBLUP
#'         $ mmsummary data summary
#' @export
#' @import dplyr
#'
#' @examples \dontrun{qcdata = dataqc(geno_data, pheno_data, qtl_data, qtl_dom_data)}
mmdata <- function(geno_data, pheno_data, qtl_data, qtl_dom_data)
{
  # Transform genotype
  rownames(geno_data) = geno_data[,2]
  mmgeno_data = t(geno_data[,-c(1:3)])

  # Summary
  traitName = colnames(pheno_data)[2]
  lineName = intersect(unique(pheno_data$GID),rownames(mmgeno_data))
  lineNum = length(lineName)

  # Prepare genotype
  mmgeno_data = as.matrix(mmgeno_data[lineName, ])

  # Prepare phenotype
  mmpheno_data = pheno_data %>%
    dplyr::filter(GID %in% lineName) %>%
    droplevels()
  colnames(mmpheno_data)[colnames(mmpheno_data) == traitName] <- deparse(substitute(trait))

  # Prepare qtl
  if(!is.null(qtl_data)){
    site_qtl = qtl_data$QTL
    m_qtl = length(site_qtl)
  } else {
    m_qtl = 0
  }

  # Prepare qtl_dom
  if(!is.null(qtl_dom_data)){
    site_qtl_dom = qtl_dom_data$QTL
    m_qtl_dom = length(site_qtl_dom)
  } else {
    m_qtl_dom = 0
  }

  # Additive relationship matrix
  A = A.mat(mmgeno_data)
  colnames(A) = rownames(A) = rownames(mmgeno_data)

  # Dominance relationship matrix
  D = D.mat(mmgeno_data)
  colnames(D) = rownames(D) = rownames(mmgeno_data)

  # Calculate reduced additive relationship matrix
  KaXa = calculateKaXa(mmgeno_data, mmpheno_data, A, site_qtl, m_qtl)

  # Calculate reduced dominance relationship matrix
  KdXd = calculateKdXd(mmgeno_data, mmpheno_data, D, site_qtl_dom, m_qtl_dom)

  return(list(mmgeno_data = mmgeno_data,
              mmpheno_data = mmpheno_data,
              Ka = KaXa$Ka,
              Xa = KaXa$Xa,
              Kd = KdXd$Kd,
              Xd = KdXd$Xd,
              A = A,
              D = D,
              mmsummary = list(traitName = traitName,
                               lineName  = lineName)))
}
```

```
#' mmgblup function
#'
#' @param data a data frame
#' @param Ka additive genetic relationship matrix
#' @param Kd dominance genetic relationship matrix
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = mmgblup(data = cbind(dt, mmdata$Xa, mmdata$Xd), Ka = mmdata$Ka, Kd = mmdata$Ka)}
mmgblup <- function(data, Ka, Kd)
{
  ## Evaluate input
  if(missing(Ka) | missing(Kd)){
    stop("Error: no input for additive kinship matrix or dominance covariance matrix")
  }

  ## Receive fixed effect column names
  ## If no fixed effect columns in data, then only intercept is fixed effect
  fixColName = names(data)[!(names(data) %in% c("GID","trait"))]
  if(length(fixColName)==0){
    fixColName = "1"
  }

  datafake = data %>% dplyr::mutate(GID1=GID)

  mod = mmer(reformulate(fixColName, "trait"),
             random = ~vsr(GID, Gu=Ka) + vsr(GID1, Gu=Kd),
             rcov = ~units,
             data = datafake,
             verbose = FALSE, date.warning = FALSE)

  BV = subset.data.frame(data, select = "GID")
  fixColName_remain = as.vector(mod$Beta$Effect)
  fixColName_remain_Aidx = grep(pattern = "_A$", x = fixColName_remain)
  fixColName_remain_Didx = grep(pattern = "_D$", x = fixColName_remain)

  BV$mu  = mod$Beta$Estimate[1]                                                                     # mu
  BV$A_l = as.matrix(data[,fixColName_remain[fixColName_remain_Aidx]]) %*%
    as.vector(mod$Beta$Estimate)[fixColName_remain_Aidx]                                            # A-major
  BV$D_l = as.matrix(data[,fixColName_remain[fixColName_remain_Didx]]) %*%
    as.vector(mod$Beta$Estimate)[fixColName_remain_Didx]                                            # D-major
  BV$A_s = mod$U$`u:GID`$trait[BV$GID]                                                              # A-minor
  BV$D_s = mod$U$`u:GID1`$trait[BV$GID]                                                             # D-minor
  BV$pre = BV$mu + BV$A_l + BV$D_l + BV$A_s + BV$D_s

  return(list(mod, BV))

}
```


```
#' pheno.generate function
#'
#' @param model a character that indicate the genetic effects included ("A" or "AD")
#' @param geno_data geno_data given from geno.generate() function
#' @param effects an M*p matrix p is the number of genetic effects included (such as "A" for 1, "AD" for 2)
#' @param sigma.error the variance for residual effects.
#'
#' @return pheno_data a data frame, with two columns $GID, and $SimTrait
#' @export
#'
#' @examples \dontrun{pheno.generate(genotypes = t(geno_data[-c(1:3)]), effects = b,
#'                    indNum = indNum, sigma.error = sigma_error)}
pheno.generate <- function(model, geno_data, effects, sigma.error){

  Ga = t(geno_data[-c(1:3)])
  indNum = nrow(Ga)

  if(model == "AD") {
    Gd = 1 - abs(Ga)
    g = tcrossprod(Ga, t(effects[,"A"])) + tcrossprod(Gd, t(effects[,"D"]))
  } else if(model == "A") {
    g = tcrossprod(Ga, t(effects[,"A"]))
  } else {
    stop("Error: Model type not recognized. Please provide either 'A' or 'AD'.")
  }

  pheno_data = data.frame(GID = as.factor(rownames(Ga)))

  g = as.vector(g)
  error = rnorm(indNum, mean = 0, sd = sqrt(sigma.error))
  pheno_data = cbind(pheno_data, g + error)

  colnames(pheno_data) = c("GID","SimTrait")

  pheno_data$GID = as.factor(pheno_data$GID)

  return(pheno_data)
}
```


```
#' qtxnetwork.output.trans
#'
#' @param pheno_data a data frame of phenotypic data, one of the outputs from simulation
#' @param pre_file path and name for .pre file
#'
#' @return a list
#'         $ qtl_data a data frame with additive qtl information
#'         $ qtl_dom_data a data frame with dominance qtl information
#' @export
#'
#' @examples \dontrun{qtl = qtxnetwork.output.trans(pheno_data, pre_file)
#'                    qtl$qtl_data
#'                    qtl$qtl_dom_data}
qtxnetwork.output.trans <- function(pheno_data, pre_file)
{
  traitName = colnames(pheno_data)[-c(1)]
  traitNum = length(traitName)

  start_title = "_1D_effect"
  end_title = "_1D_heritability"
  df.qtl = data.frame()
  for(c in 1:traitNum)
  {
    trait = traitName[c]
    file_lines = readLines(pre_file)
    # locate target lines
    start_line <- 0
    end_line <- 0
    for (i in 1:length(file_lines)) {
      line <- file_lines[i]
      if (line == start_title) {
        start_line <- i
        next
      }
      if (line == end_title) {
        end_line <- i
        break
      }
    }
    # store in data.frame
    if(start_line) {
      tmp <- as.data.frame(do.call(rbind, strsplit(file_lines[(start_line + 2):(end_line - 2)], "\\s+")))
      df.qtl = rbind(df.qtl, data.frame(trait, tmp))
    }
  }

  if(nrow(df.qtl)==0) {return(list(qtl_data = data.frame(TRAIT=character(), QTL=character()),
                                   qtl_dom_data = data.frame(TRAIT=character(), QTL=character())))}

  colnames(df.qtl) = c("TRAIT", "QTL", "SNPID", "A", "SE", "P-Value","D", "DSE", "DP-Value")

  dt = df.qtl %>% dplyr::select(c("TRAIT","SNPID","A", "D"))
  dt_reshape=reshape(dt,
                     idvar=c("TRAIT", "SNPID"),
                     varying=c("A","D"),
                     v.names="Effect",
                     timevar="Type",
                     times=c("A","D"),
                     direction="long")
  # Additive qtl
  qtl_data = dt_reshape %>%
    dplyr::filter(Type == "A") %>%
    dplyr::filter(Effect != "---") %>%
    dplyr::select(c("TRAIT", "SNPID")) %>%
    dplyr::rename(QTL = SNPID)
  # AE qtl
  qtl_dom_data = dt_reshape %>%
    dplyr::filter(Type == "D") %>%
    dplyr::filter(Effect != "---") %>%
    dplyr::select(c("TRAIT", "SNPID")) %>%
    dplyr::rename(QTL = SNPID)

  if(nrow(qtl_data)>0) rownames(qtl_data) = 1:nrow(qtl_data)
  if(nrow(qtl_dom_data)>0) rownames(qtl_dom_data) = 1:nrow(qtl_dom_data)

  return(list(qtl_data = qtl_data, qtl_dom_data = qtl_dom_data))
}
```


```
#' snp.effect function
#'
#' @param model a character that indicate the genetic effects included ("A" or "AD")
#' @param snpNum The number of SNPs
#' @param major_a_idx The index for major additive SNP
#' @param major_d_idx The index for major dominance SNP
#' @param variance_a_major The variance for major additive SNP effect
#' @param variance_a_minor The variance for minor additive SNP effect
#' @param variance_d_major The variance for major dominance SNP effect
#' @param variance_d_minor The variance for minor dominance SNP effect
#'
#' @return effect an M*p matrix p is the number of genetic effects included (such as "A" for 1, "AD" for 2)
#'
#' @export
#' @importFrom stats rnorm
#'
#' @examples snp.effect(model = "AD", snpNum = 2000,
#'                      major_a_idx = c(500, 750, 1000, 1250, 1500), variance_a_major = 0.02, variance_a_minor = 0.002,
#'                      major_d_idx = c(500, 750, 1000, 1250, 1500), variance_d_major = 0.02, variance_d_minor = 0.002)
snp.effect <- function(model, snpNum,
                       major_a_idx, variance_a_major, variance_a_minor,
                       major_d_idx, variance_d_major, variance_d_minor) {

  if(model == "A") {
    if(max(major_a_idx) > snpNum) {
      stop("Error: snp.effect(). The input major_a_idx is out of the snpNum range.")
    }

    if(missing(variance_a_major) | missing(variance_a_minor)) {
      stop("Error: snp.effect(). Please provide variance values for both major and minor additive effects.")
    }

    # determine minor index for additive effects
    minor_a_idx = setdiff(1:snpNum, major_a_idx)

    # store main effect for all SNPs
    effects = data.frame(A = rep(0, snpNum))

    # Simulate additive genetic effect
    effects[major_a_idx, "A"] <- rnorm(length(major_a_idx), mean = 0, sd = sqrt(variance_a_major))
    effects[minor_a_idx, "A"] <- rnorm(length(minor_a_idx), mean = 0, sd = sqrt(variance_a_minor))
  } else if(model == "AD") {
    if(max(major_a_idx) > snpNum | max(major_d_idx) > snpNum) {
      stop("Error: snp.effect(). The input major_a_idx or major_d_idx is out of the snpNum range.")
    }

    if(missing(variance_d_major) | missing(variance_d_minor)) {
      stop("Error: snp.effect(). Please provide variance values for both major and minor dominance effects.")
    }

    # determine minor index for additive and dominance effects
    minor_a_idx = setdiff(1:snpNum, major_a_idx)
    minor_d_idx = setdiff(1:snpNum, major_d_idx)

    # store main and interaction effect for all SNPs
    effects = data.frame(A = rep(0, snpNum), D = rep(0, snpNum))

    # Simulate genetic effects
    effects[major_a_idx, "A"] <- rnorm(length(major_a_idx), mean = 0, sd = sqrt(variance_a_major))
    effects[minor_a_idx, "A"] <- rnorm(length(minor_a_idx), mean = 0, sd = sqrt(variance_a_minor))

    effects[major_d_idx, "D"] <- rnorm(length(major_d_idx), mean = 0, sd = sqrt(variance_d_major))
    effects[minor_d_idx, "D"] <- rnorm(length(minor_d_idx), mean = 0, sd = sqrt(variance_d_minor))
  } else {
    stop("Error: snp.effect(). Unsupported model type. Please use 'A' for additive effects only or 'AD' for both additive and dominance effects.")
  }

  # Return effects matrix
  return(effects)
}
```