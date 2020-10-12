#####################################################################################################
#### -------------------------------------- main function -------------------------------------- ####
#####################################################################################################
library(progress)

MRMDR <- function(phes, snp.mat, K=2, cv=10, nperm=1000, sele.type='cvc', covrt=NULL, trim=T, test.type='ht2'){

  #adjust covariant's effect for each phenotype
  if(!is.null(covrt)) {
    fun <- function(y){
      resid <- lm(y ~ covrt)$residuals
      return(resid)
    }
    phes <- apply(phes, 2, fun)
  }

  library(fclust)
  if (trim==T){
    clust <- FKM.noise(phes, k=2, stand=1)
    data <- cbind(clust$U, phes) #***
    noise.clust <- matrix(1-clust$U[,1]-clust$U[,2],ncol=1)
    c <- cbind(clust$U, noise.clust)
    trim <- which(colnames(c)[max.col(c,ties.method="first")]=="")
    phes[trim,] <- NA
    comp <- which(complete.cases(phes))
    phes <- as.matrix(phes[comp, ])
    data <- as.matrix(data[comp, ])
    snp.mat <- as.matrix(snp.mat[comp, ])
  } else{
    clust <- FKM(phes, k=2, stand=1)
    data <- cbind(clust$U, phes) #***
    phes <- as.matrix(phes)
    data <- as.matrix(data)
    snp.mat <- as.matrix(snp.mat)
  }


  set.seed(42)
  n <- nrow(phes)
  p <- ncol(snp.mat)
  snp.combs <- combn(p, K)  ## all possible combinatory pairs
  ns <- ncol(snp.combs)
  test.stats <- rep(0L, ns)

  aa <- sample(1:n, n)  ## shuffle samples

  result <- MRMDR_cv(data[aa,], snp.mat[aa,], K, cv=10, ratio = NULL, snp.combs, sele.type ='cvc', test.type)

  model.cons <- result$cvc
  model.sele <- result$best.pair
  model.all.pair = result$all.pair # @@@
  # model.score <- result$score
  model.scores = result$scores # ***

  best.ksnps <- snp.combs[, model.sele, drop=F]
  ksnps = snp.combs[, model.all.pair, drop=F] # @@@


  # permutation test
  perm.pv <- NULL
  if(nperm > 0){

    emp_stats_null <- permutation_test(data, snp.mat[, best.ksnps], K, nperm, cv, test.type)
    # perm.pv <- mean(ifelse(emp_stats_null > model.score , 1, 0))

    n.model <- length(model.scores)
    for (i in 1:n.model){
      perm.pv[i] = mean(ifelse(emp_stats_null > model.scores[i], 1, 0))
    }
  }

  # final.result <- list(best_ksnps=best.ksnps, cvc=model.cons, score=model.score, pv= perm.pv)

  final.result <- list(best_ksnps=best.ksnps, 'ksnps'=ksnps,  # @@@
                       cvc=model.cons, scores=model.scores, pv= perm.pv)
  return(final.result)
}

#####################################################################################################
#### --------------------------------- subfunction of MRMDR -------------------------------- ####
#####################################################################################################

MRMDR_cv <- function (data, snp.all, K, cv=10, ratio = NULL, snp.combs, sele.type = 'cvc', test.type) {
  if (is.null(ratio)){
    ratio <- sum(data[, 1])/sum(data[, 2]) #***
  }

  ns <- ncol(snp.combs)
  n <- dim(data)[1]

  ## split the whole data into folds
  cvlen <- floor(n/cv)
  cc <- 1:n
  test.stats <- train.stats <- rep(0, ns)
  temptest.stats <- matrix(0, ns, cv)
  best.comb <- rep(0, cv)

  ## select best model(i.e. snp combination)
  for(i in 1:cv){
    # print(paste(i, "th cv"))
    test.ids <- ((i-1)*cvlen+1):(i*cvlen)
    train.ids <- cc[-test.ids]

    # pb_comb <- progress_bar$new(  format = "  combination [:bar] :percent eta: :eta",
    #                             total = ns, clear = FALSE, width= 60)
    # pb_comb$tick(0)
    for(j in 1:ns){
      # pb_comb$tick()
      temp.result <- MRMDR_cells(train.ids, test.ids, snp.all[, snp.combs[, j]], ratio, data, test.type)
      train.stats[j] <- temp.result$train.stat
      temptest.stats[j, i] <- temp.result$test.stat
    }

    # which snp pair has best training stat for each trainind set
    best.comb[i] <- which.max(train.stats)
  }
  test.stats <- rowMeans(temptest.stats, na.rm = TRUE)  ## average testing stats for all snp pairs

  if(sele.type == 'cvc'){
    ta <- table(best.comb)
    cvc <- ta[order(-ta)]     ## the largest cvc
    best.pair <- as.numeric(names(ta[order(-ta)]))[1]  ## the pair gets largest cvc
    all.pair = as.numeric(names(ta[order(-ta)])) # @@@
  }
  if(sele.type == 'score'){
    best.pair <- which.max(test.stats)  ## the pair gives the largest test score
    cvc <- length(which(best.comb == best.pair))
  }

  # sele.score <- c(test.stats[best.pair])
  sele.score = test.stats[all.pair] # ***

  ## sele.score -- corresponding to the test score of the final selected model
  ## test.stats -- record test.scores for all possible k-way model (snp interactions)

  ## Save the test score of the best model in each cv
  scores.cv = temptest.stats[best.pair, ]

  # return(list('cvc' = cvc, 'score' = sele.score,
  #             'best.pair' = best.pair, 'test.stats' = test.stats))

  return(list('cvc' = cvc, 'scores' = sele.score,
              'best.pair' = best.pair, 'all.pair'=all.pair, 'test.stats' = test.stats)) # ***
}


#####################################################################################################
#### ---------------------------------- subfunction of MRMDR_cv -------------------------------- ####
#####################################################################################################
MRMDR_cells <- function (train.ids, test.ids, snp.mat, ratio, data, test.type) {

  snp.mat <- as.matrix(snp.mat)

  ## remove missing values
  fids <- which(complete.cases(snp.mat))
  snp.mat <- as.matrix(snp.mat[fids, ])

  test.ids <- intersect(test.ids, fids)
  train.ids <- intersect(train.ids, fids)

  k <- ncol(snp.mat)

  ## split data into cells
  tlist <- vector('list', k)
  for(i in 1:k) tlist[[i]] <- snp.mat[, i]
  cells <- split(data.frame(cbind(fids, snp.mat)), tlist)

  ## delete NULL cells
  obs.cell <- sapply(cells, function(x) nrow(x))
  cell.null <- which(obs.cell == 0)
  if (length(cell.null) > 0){cells <- cells[-cell.null]}

  ## get trainid in each cell
  cells.trainid <- lapply(cells, function(x) return(intersect(x[, 1], train.ids)))

  cells.num <- length(cells)

  ## compare local ratio and global ratio
  high.all <- NULL
  for(i in 1:cells.num){
    temp.ids <- cells.trainid[[i]]
    if (sum(data[temp.ids, 2])==0) next #***
    if (sum(data[temp.ids, 1])/sum(data[temp.ids, 2]) >= ratio){ #***
      high.all <- c(high.all, cells[[i]][, 1])
    }
  }

  if (test.type == 'ht2'){
    train.stat <- cal_ss(train.ids, high.all, data[,-(1:2)]) #***
    test.stat <- cal_ss(test.ids, high.all, data[,-(1:2)]) #***
  }
  else if (test.type == 't'){
    train.stat = cal_tstat(train.ids, high.all, data[,-1])
    test.stat = cal_tstat(test.ids, high.all, data[,-1])
  }
  else if (test.type == 'rank_t'){
    train.stat = cal_tstat(train.ids, high.all, rank(data[,-1]))
    test.stat = cal_tstat(test.ids, high.all, rank(data[,-1]))
  }
  else if (test.type == 'mvr'){
    train.stat <- cal_mvr(train.ids, high.all, data[,-(1:2)]) #***
    test.stat <- cal_mvr(test.ids, high.all, data[,-(1:2)]) #***
  }
  else if (test.type =='uvr'){
    train.stat = cal_uvr(train.ids, high.all, data[,-1])
    test.stat = cal_uvr(test.ids, high.all, data[,-1])
  }
  else if (test.type=='rank_ht2'){
    train.stat <- cal_ss(train.ids, high.all, apply(data[,-(1:2)],2,rank))
    test.stat <- cal_ss(test.ids, high.all, apply(data[,-(1:2)],2,rank))
  }

  return(list('train.stat' = train.stat, 'test.stat' = test.stat ))
}


#####################################################################################################
#### --------------------------------- subfunction of MRMDR -------------------------------- ####
#####################################################################################################

permutation_test <- function(data, snp.all, K, nperm, cv, test.type){

  set.seed(42)
  n <- nrow(data)

  test.stats <- rep(0, cv)
  cvlen <- floor(n/cv)
  cc <- 1:n
  ratio <- sum(data[, 1])/sum(data[, 2]) #***

  run <- 0
  stats <- NULL

  pb <- progress_bar$new(  format = "  permutation test [:bar] :percent eta: :eta",
                           total = nperm, clear = FALSE, width= 60)
  pb$tick(0)
  for (i in 1:nperm){
    pb$tick()
    run <- run + 1
    perm.id <- sample(1:n, n)
    perm.phes <- data[perm.id,]

    for(j in 1:cv){
      test.ids <- ((j-1)*cvlen+1):(j*cvlen)
      train.ids <- cc[-test.ids]
      temp.result <- MRMDR_cells(train.ids, test.ids, snp.all, ratio, perm.phes, test.type)
      test.stats[j] <- temp.result$test.stat
    }

    stats <- c(stats, mean(test.stats))
  }
  return(stats)
}

#####################################################################################################
#### ------------------------------- subfunction of MRMDR_cells -------------------------------- ####
#####################################################################################################

cal_ss <- function(ids, high.all, phes){

  high.ids <- intersect(ids, high.all)
  low.ids <- setdiff(ids, high.ids)
  d <- ncol(phes)

  s1 <- phes[high.ids, ]
  s2 <- phes[low.ids, ]

  if (length(high.ids) == 0 || length(low.ids) == 0) return(0)

  s1 <- as.matrix(s1)
  s2 <- as.matrix(s2)

  if (ncol(s1) == 1){s1 <- t(s1)}
  if (ncol(s2) == 1){s2 <- t(s2)}


  stat <- dire_ht2(s1, s2, phes)$fstat                   ## another version
  ## stat is scaled that it follows a F distribution with degree d, n-1-d under the null
  return(stat)
}

## calculate HT2 directly
dire_ht2 <- function(X, Y, phes){

  # number of observations for two group:
  l1 <- nrow(X)
  l2 <- nrow(Y)
  d <- ncol(X)

  # Sample mean vectors for the each group:
  m1 <- apply(X, 2, mean)
  m2 <- apply(Y, 2, mean)

  # "pooled" sample covariance matrix:
  poolS <- ((l1-1)*cov(X)+ (l2-1)*cov(Y))/(l1+l2-2)

  if (any(is.na(poolS)) || abs(det(poolS)) < 0.00001) poolS <- cov(phes)

  # Hotelling T^2, the F-statistic, and the P-value:
  T2 <- ((l1*l2)/(l1+l2))*(t(m1-m2) %*% solve(poolS) %*% (m1-m2) )

  Fstat <- ((l1+l2-d-1)*T2)/((l1+l2-2)*d)
  # pvalue <- pf(Fstat, d, l1 + l2 - d - 1, lower.tail = FALSE)

  return(list("stat" = round(T2, 4),
              "fstat" = round(Fstat, 4) ))
}


## calculating multivariate rank score
cal_mvr <- function(ids, high.all, phes){
  # library(MNM)
  high.ids = intersect(ids, high.all)
  low.ids = setdiff(ids, high.ids)

  high.phes <- phes[high.ids, ]
  low.phes <- phes[low.ids, ]
  phes.new <- rbind(high.phes, low.phes)
  group <- factor(rep(c("high","low"), c(length(high.ids), length(low.ids))))

  if (length(high.ids) < 2 || length(low.ids) < 2) return(0)

  stat = MNM::mv.Csample.test(phes.new, group, score="r")$statistic
  return(stat)
}

## calculation univariate rank score
cal_uvr <- function(ids, high.all, phes){

  high.ids = intersect(ids, high.all)
  low.ids = setdiff(ids, high.ids)

  high.phes <- phes[high.ids]
  low.phes <- phes[low.ids]
  phes.new <- c(high.phes, low.phes)
  group <- factor(rep(c("high","low"), c(length(high.ids), length(low.ids))))

  if (length(high.ids) < 2 || length(low.ids) < 2) return(0)

  stat = kruskal.test(phes.new, group)$statistic
  return(stat)
}


## calculating t score
cal_tstat <- function(ids, high.all, phes){
  high.ids = intersect(ids, high.all)
  low.ids = setdiff(ids, high.ids)

  high.phes = phes[high.ids]
  low.phes = phes[low.ids]

  if (length(high.ids) == 0 || length(low.ids) == 0) return(0)
  if (length(union(high.ids, low.ids)) <= 2) return(0)

  stat = t.test(high.phes, low.phes, var.equal = TRUE)$statistic

  return(abs(stat))
}
