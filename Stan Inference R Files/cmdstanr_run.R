library(ggmcmc)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)
library(posterior)
library(tidyverse)
library(ggplot2)
library(plyr)
library(grid)
library(trialr)
library(ordinal)
library(cmdstanr)

options(mc.cores = parallel::detectCores())

# Theme for ggplots
theme_set(theme_bw(base_size = 16))

# Load the data
source("./modules/load_tendon_data.R")
data <- load_tendon_data(threshold = 0.3)

# Load in the tendon model
source("./modules/tendonModel_nophi_scaled_st_new.R")
stanvars <- get_model_stancode()

# Generate the model formula
source("./modules/model_formula.R")
pop_model_formula <- mixed_effects_model.nophi.st()

source("./modules/model_priors.R")
popPriors <- get_model_priors.me.nophi.rescaled.st()

# Filter the full data
type_to_infer <- c("sdft") # c("cdet")
populationData <-
  data %>% filter(type %in% type_to_infer)

args <- commandArgs(trailingOnly = T)
outputFolder <- paste0("./inference output/", args[1])
dir.create(path = outputFolder,
           recursive = T,
           showWarnings = F)

file.copy('./cmdstanr_run.R',
          paste0(outputFolder, 'script.R'),
          overwrite = T)

num_data <- length(unique(populationData$group))

# Run the Inference
num_chains <- 10
num_threads <- 1
warmup_iter <- 2e3
sampling_iter <- 4e3
num_samps <- warmup_iter + sampling_iter
num_cores <- min(num_chains * num_threads, 16)
refresh <- 1e2

stan_data <- make_standata(
  pop_model_formula,
  data = populationData,
  prior = popPriors,
  family = gaussian(),
  stanvars = stanvars
)

ids_to_remove <- c(9:21, 23, 26)
stan_data.edit <-
  stan_data[names(stan_data)[-ids_to_remove]] # Remove unused standata
names(stan_data.edit)[9:11] <- c("J", "M", "NC") # Rename standata
stan_data.edit["K"] <- num_data

source("./modules/init_vals_centred.R")
inits <-
  get_init_vals.me.nophi.rescaled.st.previous_run(num_chains, num_data)

mod <-
  cmdstan_model(paste0("./cmdstanr_stan_files/", "tendon_me.stan"),
                force_recompile = F)

fit <- mod$sample(
  data = stan_data.edit,
  init = inits,
  chains = num_chains,
  parallel_chains = num_chains,
  iter_warmup = warmup_iter,
  iter_sampling = sampling_iter,
  refresh = 1e2,
  step_size = 0.001,
  adapt_delta = 0.99,
  max_treedepth = 14
)

draws_df <-
  fit$draws(format = "df") %>% select(all_of(starts_with(params)))

groups <- unique(populationData$group)
param_names <- c("phimu", "phiE"
                 , "a", "b")

rename_group <- function(val) {
  index <- as.integer(str_extract(col_name, "[\\d]"))
  groups[index]
}

get_grouped_samples <- function(output_df) {
  get_param_samps <- function(param_name) {
    param_samps <- output_df[, grepl(param_name, colnames(output_df))]
    param_samps <- gather(param_samps, group, par, paste0(param_name, "[", 1:18, "]")) %>% select(par)
    names(param_samps)[1] <- param_name
    param_samps
  }
  
  bind_cols(lapply(param_names, get_param_samps)) %>% mutate(group = rep(groups, each =
                                                                           sampling_iter * num_chains)) %>% select(group, everything())
}

np <- nuts_params(popFit)
mcmc_nuts_energy(np) + ggtitle("NUTS Energy Diagnostic")

np %>%
  group_by(Chain) %>%
  filter(Parameter == "energy__") %>%
  summarise((sum(diff(Value) ** 2) / length(Value)) / var(Value))

divergences <- np %>% filter(Parameter == "divergent__")

summaries <-
  divergences %>% group_by(Chain) %>% summarise(num_div = sum(Value))

data <-
  data.frame(Chain = as.factor(summaries$Chain),
             num_div = summaries$num_div)

p <- ggplot(data, aes(x = Chain, y = num_div)) +
  geom_bar(stat = "identity") + ylab("Number of Divergences")
png(
  filename = paste0(outputFolder, "div_per_chain.png"),
  units = "in",
  res = 300,
  width = 6,
  height = 3
)
print(p)
dev.off()

step_sizes <- np %>% filter(Parameter == "stepsize__")
ggplot(step_sizes, aes(x=Chain, y=Value)) + geom_line()
ggsave(
  filename = paste0(outputFolder, "step_size_per_chain.png"),
  units = "in",
  dpi = 300,
  width = 6,
  height = 3
)

write.table(
  get_grouped_samples(draws_df),
  paste0(outputFolder, "output.csv"),
  sep = ",",
  row.names = F
)
saveRDS(fit, paste0(outputFolder, "fit.rds"))

fit$summary()
fit$diagnostic_summary()
fit$cmdstan_diagnose()
