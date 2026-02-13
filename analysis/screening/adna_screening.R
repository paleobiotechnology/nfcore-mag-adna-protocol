## aDNA Summary Script

## Load libraries
library(readr)
library(dplyr)


## Load data
raw <- read_tsv("data/premade_mag_results/execution-cli/bin_summary.tsv")

raw %>%
    mutate(
      is_lq_mag = if_else(Completeness_checkm < 50 & Contamination_checkm <= 10, TRUE, FALSE),
      is_mq_mag = if_else(Completeness_checkm >= 50 & Contamination_checkm <= 10, TRUE, FALSE),
        is_hq_mag = if_else(Completeness_checkm >= 90 & Contamination_checkm <= 5, TRUE, FALSE),
        is_ancient = if_else(nb_reads_aligned_pydamagebins >= 1000 & predicted_accuracy_pydamagebins >= 0.5 & qvalue_pydamagebins <= 0.05, TRUE, FALSE),

    ) |> select(bin, Completeness_checkm, Contamination_checkm, nb_reads_aligned_pydamagebins, starts_with('is_'))



