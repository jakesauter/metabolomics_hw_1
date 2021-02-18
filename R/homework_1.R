
# 1 Data preprocessing

library(tidyverse)
library(patchwork) 

data_dir <- "/home/x1/Documents/Weill_Cornell/Functional_Interpretation/metabolomics/homework_1/data/"
data_filename <- list.files(data_dir, "Simulated_metabolomics_data.RData", 
                            full.names = TRUE)
data_filename %<>% normalizePath(mustWork=TRUE)
data_env <- new.env()
load(data_filename, envir = data_env)

## Data now loaded, checking for data quality 

plot_list <- vector('list', 9)

df <- 
  data_env %>% 
  .$dat %>%
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-c(sample_num), 
               names_to='metabolite', values_to="concentration")

for (group_num in 1:9) {
  first_sample <- 1 + (100 * (group_num - 1))
  last_sample <- 100 * (group_num - 1) + 100
  
  plot_list[[group_num]] <- df %>% 
    filter(sample_num >= first_sample, 
           sample_num <= last_sample) %>% 
    ggplot() + 
      geom_boxplot(aes(x=sample_num, 
                       y=concentration, 
                       group=sample_num)) + 
      ylim(c(-4, 6))
}


(plot_list[[1]] + plot_list[[2]] + plot_list[[3]]) / 
(plot_list[[4]] + plot_list[[5]] + plot_list[[6]]) / 
(plot_list[[7]] + plot_list[[8]] + plot_list[[9]])

# we find that some samples have much higher standard
# measured concentration levels, lets try to normalize


logged_data <- log(data_env$dat)
mean(logged_data)
sd(logged_data)

data_env$dat %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-sample_num) %>% 
  ggplot() + 
  geom_histogram(aes(x=value), bins=30)


#################################
# Quotient Normalization 
#################################
# 1. find the median sample
# 2. calculate the quotient of every
#    sample with the median sample 
# 3. Find the median of these quotients
# 4. divide all variables by this median
#################################

median_sample <- apply(data_env$dat, 2, median)
normalized_dat <- data_env$dat 

for (i in 1:nrow(normalized_dat)) {
  quotient <- data_env$dat[i,] / median_sample
  median_quotient <- median(quotient)
  normalized_dat[i,] = data_env$dat[i,] / median_quotient
}


##############################################



df_before <- 
  data_env %>% 
  .$dat %>%
  .[1:100, ] %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-c(sample_num), 
               names_to='metabolite', values_to="concentration")


df_after <- 
  normalized_dat %>%
  .[1:100, ] %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-c(sample_num), 
               names_to='metabolite', values_to="concentration")

p1 <- 
  df_before %>% 
  ggplot() + 
  geom_boxplot(aes(x=sample_num, 
                   y=concentration, 
                   group=sample_num))

p2 <- 
  df_after %>% 
  ggplot() + 
  geom_boxplot(aes(x=sample_num, 
                   y=concentration, 
                   group=sample_num))

p1 / p2


# 2 Phenotype association

?t.test

p_vals <- vector('double', ncol(dat))
for (metabolite in 1:ncol(dat)) {
  print(metabolite)
  male <- dat[gender==1, metabolite]
  female <- dat[gender==2, metabolite]
  p_vals[metabolite] <- t.test(male, female)$p.value
}

df <- data.frame(
        metabolite=colnames(data_env$dat), 
        pval=p_vals) 

df %>% 
  ggplot() + 
  geom_histogram(aes(x=pval))

df %>% 
  filter(p_vals < 0.05) %>% 
  ggplot() + 
    geom_histogram(aes(x=pval)) 

df %>% 
  filter(p_vals < 0.05) %>% 
  nrow()


# Plot p-values as a distribution (use plot of your choice).
# What do you observe? Perform multiple testing correction 
# with a method and α of your choice. How many metabolites
# are significantly associated with gender before and after
# multiple testing correction?

hist(p_vals, 30)

data_env$dat %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-sample_num) %>% 
  ggplot() + 
  geom_histogram(aes(x=value), bins=30)

normalized_dat %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-sample_num) %>% 
  ggplot() + 
  geom_histogram(aes(x=value), bins=30)


# Which 10 metabolites are most significantly associated with gender? 
# Which pathways do these metabolites belong to? Visualize the 
# differences between male and female for the top hit.

# 3 Gaussian graphical model (GGM)

plot_list <- vector('list', 9)

df <- 
  data_env %>% 
  .$dat %>%
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-c(sample_num), 
               names_to='metabolite', values_to="concentration")

for (group_num in 1:9) {
  first_sample <- 1 + (100 * (group_num - 1))
  last_sample <- 100 * (group_num - 1) + 100
  
  plot_list[[group_num]] <- df %>% 
    filter(sample_num >= first_sample, 
           sample_num <= last_sample) %>% 
    ggplot() + 
    geom_boxplot(aes(x=sample_num, 
                     y=concentration, 
                     group=sample_num)) + 
    ylim(c(-4, 6))
}


(plot_list[[1]] + plot_list[[2]] + plot_list[[3]]) / 
  (plot_list[[4]] + plot_list[[5]] + plot_list[[6]]) / 
  (plot_list[[7]] + plot_list[[8]] + plot_list[[9]])




logged_data <- log(data_env$dat)
mean(logged_data)
sd(logged_data)


data_env$dat %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-sample_num) %>% 
  ggplot() + 
  geom_histogram(aes(x=value), bins=30)


# We will first determine the **median of all samples**
# to use as our reference sample during Quotient Normalization. 

median_sample <- apply(data_env$dat, 2, median)

normalized_dat <- data_env$dat

for (i in 1:nrow(normalized_dat)) {
  normalized_dat[i,] = median_sample / data_env$dat[i,] 
}

logged_data <- log(normalized_dat)
mean(logged_data)
sd(logged_data)



# TODO: Not sure about this step yet: 
# **normalized_dat <- normalized_dat / mean(normalized_dat)**



normalized_dat %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-sample_num) %>% 
  ggplot() + 
  geom_histogram(aes(x=value), bins=30)

library(tidyverse)
library(patchwork)

plot_list <- vector('list', 9)

df_before <- 
  data_env %>% 
  .$dat %>%
  .[1:30, ] %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-c(sample_num), 
               names_to='metabolite', values_to="concentration")


df_after <- 
  normalized_dat %>%
  .[1:30, ] %>% 
  log() %>% 
  as.data.frame() %>% 
  mutate(sample_num = row_number()) %>% 
  pivot_longer(-c(sample_num), 
               names_to='metabolite', values_to="concentration")

p1 <- 
  df_before %>% 
  ggplot() + 
        geom_boxplot(aes(x=sample_num, 
                         y=concentration, 
                         group=sample_num))

p2 <- 
  df_after %>% 
  ggplot() + 
        geom_boxplot(aes(x=sample_num, 
                         y=concentration, 
                         group=sample_num))

p1 / p2


## **Phenotype association**

### **Male vs Female **


metabolites <- colnames(normalized_dat)
p_vals <- vector('double', length(metabolites))
names(p_vals) <- metabolites
for (metabolite in metabolites) {
  male   <- normalized_dat[data_env$gender==1, metabolite]
  female <- normalized_dat[data_env$gender==2, metabolite]
  p_vals[metabolite] <- t.test(male, female)$p.value
}

df <- 
  data.frame(
    metabolite=colnames(data_env$dat), 
    pval=p_vals
  ) 

df %>% 
  ggplot() + 
  geom_histogram(aes(x=pval))

df %>% 
  filter(p_vals < 0.05) %>% 
  ggplot() + 
  geom_histogram(aes(x=pval))




########################################
data_dir <- "/home/x1/Documents/Weill_Cornell/Functional_Interpretation/metabolomics/homework_1/data/"
data_filename <- list.files(data_dir, "Simulated_metabolomics_data.RData", 
                            full.names = TRUE)
data_filename %<>% normalizePath(mustWork=TRUE)
data_env <- new.env()
load(data_filename, envir = data_env)

median_sample <- apply(data_env$dat, 2, median)
normalized_dat <- data_env$dat 

for (i in 1:nrow(normalized_dat)) {
  quotient <- data_env$dat[i,] / median_sample
  median_quotient <- median(quotient)
  normalized_dat[i,] = data_env$dat[i,] / median_quotient
}

log_norm_dat <- log(normalized_dat)

metabolites <- colnames(log_norm_dat)
p_vals <- vector('double', length(metabolites))
names(p_vals) <- metabolites
for (metabolite in metabolites) {
  male   <- log_norm_dat[data_env$gender==1, metabolite]
  female <- log_norm_dat[data_env$gender==2, metabolite]
  p_vals[metabolite] <- t.test(male, female)$p.value
}

## Tell R that we would like to plot
# 1 row and 2 columns
par(mfrow = c(1, 2))

sig_p_vals <- p_vals[p_vals < 0.05]
hist(sig_p_vals, 
     main='Significant Original P-values', 
     xlab='P values')

adj_p_vals <- p.adjust(p_vals, 'bonferroni')
sig_adj_p_vals <- adj_p_vals[adj_p_vals < 0.05]
hist(sig_adj_p_vals, 
     main='Significant Bonferroni corrected P-values', 
     xlab='P values')

adj_p_vals <- p.adjust(p_vals, 'bonferroni')
length(adj_p_vals[adj_p_vals] < 0.05)

#################################################################3

before <- length(which(p_vals < 0.05))
after <- length(which(adj_p_vals < 0.05))
barplot(c(before, after), 
        main = "Significant Metabolites by Sex", 
        ylab = "Number of Significant Metabolites")


##########################################33

top_hit <- names(sort(adj_p_vals)[1])

male <- log_norm_dat[data_env$gender==1, top_hit]
female <- log_norm_dat[data_env$gender==2, top_hit]


hist(male, 
     main = paste('Male Concentrations of', top_hit), 
     col=rgb(1,0,0,0.5), 
     xlab = "Male", 
     xlim=c(-4,4), 
     ylim=c(0,150))

hist(female, 
     main = paste('Female Concentrations of', top_hit), 
     col = rgb(0,0,1,0.5), 
     xlab = "Female", 
     add=TRUE)

#############################################3

top_hit <- names(sort(adj_p_vals)[1])

male <- log_norm_dat[data_env$gender==1, top_hit]
female <- log_norm_dat[data_env$gender==2, top_hit]

male <- c(male, rep(NA, length(female) - length(male)))

df <- 
  data.frame(Male=male, 
             Female=female) %>%
  mutate(index = row_number()) %>% 
  pivot_longer(-c(index), 
               names_to='sex', values_to="concentration") %>% 
  mutate(sex = factor(sex, levels = c('Male', 'Female')))
  


df %>% 
  ggplot() + 
  geom_boxplot(aes(x=sex, 
                   y=concentration, 
                   fill=sex),
               ) + 
  ylim(c(-4, 6)) + 
  theme(legend.position = "none") + 
  labs(title=paste('Male vs Female Concentration of', top_hit)) + 
  xlab("") + 
  ylab("Normalized Concentration")
  

###################################

data_dir <- "/home/x1/Documents/Weill_Cornell/Functional_Interpretation/metabolomics/homework_1/data/"
data_filename <- list.files(data_dir, "Simulated_metabolomics_data.RData", 
                            full.names = TRUE)
data_filename %<>% normalizePath(mustWork=TRUE)
data_env <- new.env()
load(data_filename, envir = data_env)

median_sample <- apply(data_env$dat, 2, median)
normalized_dat <- data_env$dat 

for (i in 1:nrow(normalized_dat)) {
  quotient <- data_env$dat[i,] / median_sample
  median_quotient <- median(quotient)
  normalized_dat[i,] = data_env$dat[i,] / median_quotient
}

log_norm_dat <- log(normalized_dat)

n <- ncol(log_norm_dat)
n_tests <- (n*(n-1))/2 
pcor <- ppcor::pcor(log_norm_dat)
p_vals <- pcor$p.value
adj_p_vals <- p_vals * n_tests
adj_p_vals[adj_p_vals > 1] = 1
hist(adj_p_vals[adj_p_vals < 0.05])

adj_mat <- pcor$estimate
adj_mat[adj_p_vals > 0.01] = 0

for (i in 1:nrow(adj_mat)) {
  adj_mat[i, i] = 0
}

library(igraph)
graph <- graph.adjacency(adj_mat, mode = 'undirected', diag = FALSE)
x11()
plot(graph)









