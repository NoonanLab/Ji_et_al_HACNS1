library(readr)
library(dplyr)
library(tibble)
library(igraph)
library(purrr)
library(tidyr)

FL_background <- read_csv("palmer_scratch/network_deg/try/background/FL_background.csv")
HL_background <- read_csv("palmer_scratch/network_deg/try/background/HL_background.csv")
PA1_background <- read_csv("palmer_scratch/network_deg/try/background/PA1_background.csv")
PA2_background <- read_csv("palmer_scratch/network_deg/try/background/PA2_background.csv")

FL_DEG <- read_csv("~/palmer_scratch/network_deg/try/FL_DEG.csv")
HL_DEG <- read_csv("~/palmer_scratch/network_deg/try/HL_DEG.csv")
PA1_DEG <- read_csv("~/palmer_scratch/network_deg/try/PA1_DEG.csv")
PA2_DEG <- read_csv("~/palmer_scratch/network_deg/try/PA2_DEG.csv")

FL_net <- read_csv("~/palmer_scratch/network_deg/try/FL_network.csv")
HL_net <- read_csv("~/palmer_scratch/network_deg/try/HL_network.csv")
PA1_net <- read_csv("~/palmer_scratch/network_deg/try/PA1_network.csv")
PA2_net <- read_csv("~/palmer_scratch/network_deg/try/PA2_network.csv")

PA1_Gbx2_direct <- PA1_net %>%
  filter(TF == "Gbx2") %>%
  pull(target) %>%
  unique()

PA2_Gbx2_direct <- PA2_net %>%
  filter(TF == "Gbx2") %>%
  pull(target) %>%
  unique()

FL_Gbx2_direct <- FL_net %>%
  filter(TF == "Gbx2") %>%
  pull(target) %>%
  unique()

HL_Gbx2_direct <- HL_net %>%
  filter(TF == "Gbx2") %>%
  pull(target) %>%
  unique()

PA_shared_Gbx2_direct <- intersect(PA1_Gbx2_direct, PA2_Gbx2_direct)
Limb_shared_Gbx2_direct <- intersect(FL_Gbx2_direct, HL_Gbx2_direct)

high_confidence_gbx2_direct_targets <- function(net, shared_direct_targets) {
  net %>%
    filter(
      TF != "Gbx2" |
        (TF == "Gbx2" & target %in% shared_direct_targets)
    )
}

PA1_net_hc <- high_confidence_gbx2_direct_targets(PA1_net, PA_shared_Gbx2_direct)
PA2_net_hc <- high_confidence_gbx2_direct_targets(PA2_net, PA_shared_Gbx2_direct)

FL_net_hc <- high_confidence_gbx2_direct_targets(FL_net, Limb_shared_Gbx2_direct)
HL_net_hc <- high_confidence_gbx2_direct_targets(HL_net, Limb_shared_Gbx2_direct)

count_distance <- function(network_df, background_genes, source_gene = "Gbx2") {
  g <- graph_from_data_frame(network_df %>% select(TF, target), directed = TRUE)
  
  out <- tibble(
    gene = background_genes,
    distance_from_Gbx2 = NA,
    step_group = "not within 3 steps"
  )
  
  dist_vec <- distances(g, v = source_gene, to = V(g), mode = "out")[1, ]
  
  dist_df <- tibble(gene = names(dist_vec), distance_from_Gbx2 = as.numeric(dist_vec)) %>% filter(is.finite(distance_from_Gbx2)) %>% filter(distance_from_Gbx2 > 0)
  
  out <- out %>%
    select(gene) %>%
    left_join(dist_df, by = "gene") %>%
    mutate(
      step_group = case_when(
        distance_from_Gbx2 == 1 ~ "1 step",
        distance_from_Gbx2 == 2 ~ "2 steps",
        distance_from_Gbx2 == 3 ~ "3 steps",
        TRUE ~ "not within 3 steps"
      )
    )
  
  out
}

PA1_distance <- count_distance(PA1_net_hc, PA1_background)
PA2_distance <- count_distance(PA2_net_hc, PA2_background)
FL_distance  <- count_distance(FL_net_hc,  FL_background)
HL_distance  <- count_distance(HL_net_hc,  HL_background)

run_permutation <- function(distance, deg, tissue, n_perm = 20000, seed = 123) {
  set.seed(seed)
  
  shifted_genes <- deg %>% pull(1) %>% as.character()
  
  step_levels <- c("1 step", "2 steps", "3 steps", "not within 3 steps")
  
  bg <- distance %>%mutate(step_group = factor(step_group, levels = step_levels))
  
  background_genes <- unique(bg$gene)
  shifted_genes <- intersect(shifted_genes, background_genes)
  n_shifted <- length(shifted_genes)
  
  observed <- bg %>%
    mutate(shifted = gene %in% shifted_genes) %>%
    group_by(step_group) %>%
    summarise(
      observed_n = sum(shifted),
      total_background_n = n(),
      .groups = "drop"
    ) %>%
    complete(
      step_group = factor(step_levels, levels = step_levels),
      fill = list(observed_n = 0, total_background_n = 0)
    )
  
  gene_to_step <- bg$step_group
  names(gene_to_step) <- bg$gene
  
  perm_mat <- replicate(n_perm, {
    fake_shifted <- sample(background_genes, size = n_shifted, replace = FALSE)
    tab <- table(factor(gene_to_step[fake_shifted], levels = step_levels))
    as.numeric(tab)
  })
  
  perm <- as.data.frame(t(perm_mat))
  colnames(perm) <- step_levels
  
  result <- observed %>%
    mutate(
      tissue = tissue,
      n_shifted = n_shifted
    ) %>%
    rowwise() %>%
    mutate(
      expected_mean = mean(perm[[as.character(step_group)]]),
      expected_sd = sd(perm[[as.character(step_group)]]),
      fold_over_expected = observed_n / expected_mean,
      
      p_upper = sum(perm[[as.character(step_group)]] >= observed_n) / (n_perm),
      p_lower = sum(perm[[as.character(step_group)]] <= observed_n) / (n_perm),
      
      empirical_p = min(1, 2 * min(p_upper, p_lower))
    ) %>%
    ungroup() %>%
    mutate(
      p_adj = p.adjust(empirical_p, method = "BH"),
      direction = case_when(
        fold_over_expected > 1 ~ "enriched",
        fold_over_expected < 1 ~ "depleted",
        TRUE ~ "as_expected"
      )
    )
  
  list(
    observed_expected = result,
    permutations = perm
  )
}

PA1_perm <- run_permutation(
  distance = PA1_distance,
  deg = PA1_DEG,
  tissue = "PA1",
  n_perm = 20000,
  seed = 101
)

PA2_perm <- run_permutation(
  distance_df = PA2_distance,
  deg = PA2_DEG,
  tissue = "PA2",
  n_perm = 20000,
  seed = 102
)

FL_perm <- run_permutation(
  distance_df = FL_distance,
  deg = FL_DEG,
  tissue = "FL",
  n_perm = 20000,
  seed = 103
)

HL_perm <- run_permutation(
  distance_df = HL_distance,
  deg = HL_DEG,
  tissue = "HL",
  n_perm = 20000,
  seed = 104
)