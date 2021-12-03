# Code from Steve Rozen, 2021 01 24
# Analyses for Wu Yang's paper, "Evaluation of Mutational
# Signature Software on Correlated Signatures, Wu et al."
#
# This code assembles run-level results that are in individual .csv files
# into single tibbles - all_results.csv,
# and then does two levels of summary:
#
#  - at the level of approach, data set and runs, and
# stats_by_approach_and_data_set.csv - at the level of approach and data set (Supp Table S5 and S8), and
# stats_by_approach.csv - at the level of approach (Supp Table S)

library(tibble)
library(data.table)
library(dplyr)
library("readr")

myread <- function(ff) as_tibble(fread(ff)[ , -1])

# These are the files to be logically "cbind"ed

fff <- c("compositeMeasure.csv",
         "cossim.to.SBS1.csv",
         "cossim.to.SBS5.csv",
         "PPV.csv",
         "TPR.csv",
         "num.sigs.similar.to.SBS1.csv",
         "num.sigs.similar.to.SBS5.csv",
         "NumSigsExtracted.csv",
         "falseNeg.csv",
         "falsePos.csv")

dd <- lapply(fff, myread)

# Check to make sure the rows in all the files were
# for the same approach, dataset, and run
com <- dd[[1]]
for (ii in c(1, 3:5)) {
  cat(
    all(
      unlist(
        lapply(dd, function(xx, col) all(com[ ,  col] == xx[ , col]), ii))),
    "\n")
}

# Make one big tibble (cc)
cc <- dd[[1]]
for (ii in 2:length(dd)) {
  cc2 <- full_join(cc, dd[[ii]])
  cc <- cc2
}

# Make the column names more user-friendly
ccc <- cc[ , c(3, 5, 4, 1, 2, 6:14)]
colnames(ccc) <- c("Approach",
                   "SBS1_SBS5_ratio",
                   "Rsq",
                   "Run",
                   "Composite",
                   "Sim_to_SBS1",
                   "Sim_to_SBS5",
                   "PPV",
                   "TPR",
                   "N_SBS1_like",
                   "N_SBS5_like",
                   "N_sigs",
                   "FN",
                   "FP")

# 1. Write the detail table
readr::write_csv(ccc, "all_results.csv")

# 2. Summarize by approach and data set, and
# re-order / rename the columns
by.approach.and.data.set <- ccc %>%
  group_by(Approach, SBS1_SBS5_ratio, Rsq) %>%
  summarise_at(.vars = colnames(ccc)[5:14], .funs = c(mean, sd))

# Name-fixing function for summary tibbles with two-statistics (mean, sd)
fix.names <- function(nn) {
  vv <- sub("(.*)_fn1", "mean(\\1)", x = nn, perl = TRUE)
  sub("(.*)_fn2", "sd(\\1)", x = vv, perl = TRUE)
}

# Name-fixing function for summary tibbles with three-statistics (mean, median sd)
fix.names.3stats <- function(nn) {
  vv <- sub("(.*)_fn1", "mean(\\1)", x = nn, perl = TRUE)
  ww <- sub("(.*)_fn2", "median(\\1)", x = vv, perl = TRUE)
  sub("(.*)_fn3", "sd(\\1)", x = ww, perl = TRUE)
}


colnames(by.approach.and.data.set) <- fix.names(colnames(by.approach.and.data.set))

perm <- c(1:3, unlist(lapply(4:13, function(x) c(x, x + 10))))
colnames(by.approach.and.data.set)[perm] # check
by.approach.and.data.set <- by.approach.and.data.set[ , perm]

readr::write_csv(by.approach.and.data.set, "stats_by_approach_and_data_set.csv")

# 3. Summarize by approach
#
# This one is for checking supplementary tables, which were assembled
# separately.  Summarize by approach, and re-order / rename the
# columns
by.approach <- ccc %>%
  group_by(Approach) %>%
  summarise_at(.vars = colnames(ccc)[5:14], .funs = c(mean, median, sd))
colnames(by.approach) <- fix.names.3stats(colnames(by.approach))


perm3 <- c(1, unlist(lapply(2:11, function(x) c(x, x+10, x+20))))
colnames(by.approach)[perm3] # check
by.approach <- by.approach[ , perm3]

readr::write_csv(by.approach, "stats_by_approach.csv")

# 4. Wilcoxon rank sum test between approaches
#
# This test is to evaluate whether the composite measure between approaches
# are significantly different.

approaches <- ccc %>% select("Approach") %>% unique() %>% unlist() %>% unname()
index <- length(approaches)

pairwise.signif <- data.frame(
  approach1 = character(0),
  approach2 = character(0),
  p.value = numeric(0))


for(ii in seq(1,index-1)){

  for (jj in seq(ii+1,index)){

    comp <- ccc %>%
      select(Approach,Composite)

    compII <- comp %>%
      filter(Approach %in% approaches[ii]) %>%
      select(Composite) %>%
      unlist() %>% unname()

    compJJ <- comp %>%
      filter(Approach %in% approaches[jj]) %>%
      select(Composite) %>%
      unlist() %>% unname()

    res <- stats::wilcox.test(compII,compJJ)

    current <- data.frame(
      approach1 = approaches[ii],
      approach2 = approaches[jj],
      p.value = res$p.value)

    pairwise.signif <- rbind(pairwise.signif,current)

  }

}

readr::write_csv(pairwise.signif,"Pairwise.Wilcoxon.Rank.Sum.csv")
