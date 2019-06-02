test_BDMMA=function(){
  counts <- rmultinom(20,100,rep(0.2,5))
  main <- rbinom(20,1,0.5)
  confounder <- rnorm(20,0,1)
  batch <- c(rep(1,10),rep(2,10))

  library(SummarizedExperiment)
  col_data <- DataFrame(main, confounder, batch)
  mcols(col_data)$continous <- c(0L, 1L, 0L)
  ### pack different datasets into a SummarizedExperiment object
  Microbiome_dat <- SummarizedExperiment(list(counts), colData=col_data)

  output <- BDMMA(Microbiome_dat = Microbiome_dat, burn_in = 100, sample_period = 100)
}
