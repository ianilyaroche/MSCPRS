rm(list = ls()); library(tidyverse); library(car); library(ggrepel); library(sandwich); library(lmtest)

# Inputs & Outputs Folder
trial_data <- read.csv("C:/Users/L042952/Desktop/Clean_Data.csv",header = TRUE, stringsAsFactors = FALSE)
output_dir <- "C:/Users/L042952/Desktop"

# For categorical data, for modelling, requires creating factors, and uses a reference
trial_data$Phases <- relevel(factor(trial_data$Phases), ref = "PHASE3")
trial_data$disease_category <- relevel(factor(trial_data$disease_category), ref = "Cardiometabolic")
trial_data$intervention_type <- relevel(factor(trial_data$intervention_type),ref = "DRUG")

# Country inclusion threshold
region_columns <- c("has_NA", "has_LATAM", "has_EMEA", "has_ASIAPAC")
all_country_columns <- setdiff(grep("^has_", names(trial_data), value = TRUE),region_columns)
trial_counts <- sort(sapply(all_country_columns, function(country_column) sum(trial_data[[country_column]])), decreasing = TRUE)
included_country_columns <- names(trial_counts[trial_counts >= 1000])
length(included_country_columns)

# Country count below the 1,000 trial treshold
below_threshold_columns <- setdiff(all_country_columns,included_country_columns)
length(below_threshold_columns)

# Summary count of those above 1K trials
trial_data$has_Other_country <- as.integer(rowSums(trial_data[, below_threshold_columns, drop = FALSE]) > 0)
print(trial_counts[trial_counts >= 1000])

# Master function forumla for cluster robust estimates 
cluster_robust_estimates <- function(model, data, column_name = "Sponsor", confidence_level = 0.95) 
{
  
# Identify the same Sponsors for each trial used in the model
model_rows <- rownames(model.frame(model))
data_rows <- match(model_rows, rownames(data))
cluster_ids <- data[[column_name]][data_rows]
  
# Number of unique clusters of Sponsors
n_clusters <- length(unique(cluster_ids))
  
# Cluster robust regression results (ref ?coeftest and ?vcovCL)
clustered_results <- coeftest(model, vcov. = vcovCL(model,cluster = cluster_ids, type = "HC1"), df = n_clusters - 1)
  
# The critical value of t for the confidence interval
t_critical_value <- qt(1 - (1 - confidence_level) / 2, df = n_clusters - 1)
  
# Table output
tibble(term = rownames(clustered_results),
estimate = clustered_results[, 1],
std.error = clustered_results[, 2],
statistic = clustered_results[, 3],
p.value = clustered_results[, 4],
conf.low = clustered_results[, 1] - t_critical_value * clustered_results[, 2],
conf.high = clustered_results[, 1] + t_critical_value * clustered_results[, 2],
n_clusters = n_clusters)
}

# Quick references
rq1_countries <- setdiff(included_country_columns,"has_United_States")
country_string <- paste(rq1_countries, collapse = " + ")

# RQ1 model Equation
rq1_formula <- as.formula(paste( "log_duration ~ Phases + log_enrollment + disease_category +",
"intervention_type + n_sites + is_multiregional + start_year +", country_string,"+ has_Other_country"))

rq1_model <- lm(rq1_formula,data = trial_data)
print(summary(rq1_model))

# RQ1 with clusterfunction applies
rq1_model_sponsorclusters <- cluster_robust_estimates(rq1_model, trial_data)

# Clean up on labelling from columns for plot labelling
format_covariate_label <- function(term_names) {
is_intervention <- grepl("^intervention_type", term_names)
label <- gsub("^Phases|^disease_category|^intervention_type","", term_names)
label <- gsub("log_enrollment", "Log enrolment", label)
label <- gsub("^n_sites$", "Number of sites (per 10)", label)
label <- gsub("is_multiregional", "Multiregional conduct", label)
label <- gsub("start_year","Start year", label)
label <- gsub("_"," ", label)
label <- gsub("PHASE ?", "Phase ",label)
label <- gsub("EARLY","Early",label)
label <- gsub("\\|","/",label)
label <- trimws(label)

# Intervention category uses caps by default from the raw export, so this normalises text, for example 'DRUG' to 'Drug' 
label[is_intervention] <- paste0(substring(label[is_intervention], 1, 1), tolower(substring(label[is_intervention], 2)))
  
# There is 'Other' for disease_category  and intervention_type, thus a distinguish label
label[is_intervention & label == "Other"] <- "Other (intervention)" 
return(label)
}

# RQ1 model contains countries, removed via selection as plot is for characteristics only, then just a per site from 1 to 10 (rescale) - rest is 1
covariates <- rq1_model_sponsorclusters %>% filter(!grepl("^has_", term), term != "(Intercept)") %>% mutate(scale_factor = ifelse(term == "n_sites", 10, 1),

# Calculation for expressing CIs and Coef in percentage (not month)                                                                                                            
percent_change = (exp(estimate * scale_factor) - 1) * 100,
percent_change_low = (exp(conf.low * scale_factor) - 1) * 100,
percent_change_high = (exp(conf.high * scale_factor) - 1) * 100)

# Table output
covariate_data <- covariates %>% mutate(label = format_covariate_label(term),
tone = factor(ifelse( p.value >= 0.05, "Not significant", "Significant (p < 0.05)"),
levels = c("Significant (p < 0.05)", "Not significant"))) %>% arrange(percent_change) %>%
mutate(label = factor(label, levels = unique(label)))
print(covariates %>% mutate(label = format_covariate_label(term)) %>%
select(label, percent_change, percent_change_low, percent_change_high, p.value), n = 40)

# Plotting fig 4
figure4 <- ggplot(covariate_data, aes(percent_change, label, colour = tone)) + theme_minimal(base_size = 11) +
geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
geom_errorbar(aes(xmin = percent_change_low, xmax = percent_change_high),
orientation = "y", width = 0, linewidth = 0.7, alpha = 0.6) + geom_point(size = 2.6) +
scale_colour_manual(values = setNames(c("#2C7BB6", "#9E9E9E"),c("Significant (p < 0.05)", "Not significant")),
drop = FALSE, name = NULL) +

# Scaling for both negative and positive percentage changes
scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 12, base = exp(1)),
breaks = c(-40,-20,-10,0,10,20,50,100,200),
labels = function(x) {paste0(ifelse(x > 0, "+", ""), x, "%")},
expand = expansion(mult = c(0.02, 0.08))) +
labs(x = "shorter  ·  % change in trial duration vs reference (95% CI)  ·  longer", y = NULL,
title = "Structural Predictors of Trial Duration",
subtitle = paste0("References: Phase 3, Drug, Single region, Cardiometabolic area", "  | Adjusted R² = ",
round(summary(rq1_model)$adj.r.squared,3))) +
  
theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
panel.grid.minor.x = element_blank(), plot.title = element_text(face = "bold"),
plot.subtitle = element_text(colour = "grey35",size = 9),
plot.margin = margin(t = 5.5, r = 12, b = 5.5,l = 5.5) )

print(figure4)
ggsave(file.path(output_dir,"Figure4 Trial Characteristics.png" ),
figure4, width = 9, height = 10.5, dpi = 300)
cat("Complete")