rm(list = ls()); library(tidyverse); library(car); library(ggrepel); library(sandwich); library(lmtest)

#Input & Output
trial_data <- read.csv( "C:/Users/L042952/Desktop/Clean_Data.csv", header = TRUE, stringsAsFactors = FALSE)
output_dir <- "C:/Users/L042952/Desktop"

# Converts categories into levels so they can work in the model maths
trial_data$Phases <- relevel(factor(trial_data$Phases), ref = "PHASE3")
trial_data$disease_category <- relevel(factor(trial_data$disease_category), ref = "Cardiometabolic")
trial_data$intervention_type <- relevel(factor(trial_data$intervention_type), ref = "DRUG")

# Function for coefficien calculation formula for plots & tables 
to_months <- function(log_coef, baseline) {baseline * (exp(log_coef) - 1) }
to_percent <- function(log_coef) {(exp(log_coef) - 1) * 100}

# Country inclusion threshold (Removes the region columns and leaves only the countries, sets 1,000 tresh)
region_columns <- c("has_NA", "has_LATAM", "has_EMEA", "has_ASIAPAC")
all_country_columns <- setdiff(grep("^has_", names(trial_data), value = TRUE), region_columns)
trial_counts <- sort(sapply(all_country_columns, function(country_column) sum(trial_data[[country_column]])), decreasing = TRUE)
included_country_columns <- names(trial_counts[trial_counts >= 1000]) ; length(included_country_columns)
below_threshold_columns <- setdiff(all_country_columns, included_country_columns); length(below_threshold_columns)
trial_data$has_Other_country <- as.integer(rowSums(trial_data[, below_threshold_columns, drop = FALSE]) > 0)
print(trial_counts[trial_counts >= 1000])

# Quick references
rq1_countries <- setdiff(included_country_columns, "has_United_States") # Removes US country from model
country_string <- paste(rq1_countries, collapse = " + ")
n_country_tests <- length(rq1_countries) # for bonferroni factor
us_baseline <- median(trial_data$duration_months[trial_data$has_United_States == 1]); round(us_baseline,1)
germany_baseline <- median(trial_data$duration_months[trial_data$has_Germany == 1]); round(germany_baseline, 1)
cat("Sponsors", length(unique(trial_data$Sponsor)), "| Largest sponsor contributes", max(table(trial_data$Sponsor)), "trials")

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
n_clusters = n_clusters)}

# RQ1 model Equation
rq1_formula <- as.formula(paste("log_duration ~ Phases + log_enrollment + disease_category + intervention_type +", "n_sites + is_multiregional + start_year +", country_string, "+ has_Other_country"))
rq1_model <- lm(rq1_formula, data = trial_data)
rq1_model_summary   <- summary(rq1_model); print(rq1_model_summary)

# RQ1 model with the cluster robust estimates 
rq1_model_sponsorclusters  <- cluster_robust_estimates(rq1_model, trial_data)

# Table of outputs
countries <- rq1_model_sponsorclusters %>%
filter(grepl("^has_", term)) %>%
mutate(country     = gsub("_", " ", gsub("has_", "", term)), # will make has_Other_Country as 'Other Country'
months = to_months(estimate,  us_baseline),
ci_low = to_months(conf.low,  us_baseline),
ci_high = to_months(conf.high, us_baseline),
is_country = country != "Other country", # Puts Other country as false - excluded from calculation
n_trials = sapply(term, function(each_term) sum(trial_data[[each_term]])), # count of trials per country indicator
p_bonferroni = ifelse(is_country, pmin(p.value * n_country_tests, 1), NA_real_),
tier = factor(ifelse(!is_country | p_bonferroni >= 0.05, "Not significant", "Significant after Bonferroni correction"),
levels = c("Significant after Bonferroni correction", "Not significant"))) %>% arrange(months)

print(countries %>% transmute(country, n_trials, months = round(months, 2),
ci_low = round(ci_low, 2), ci_high = round(ci_high, 2),p_uncorrected = signif(p.value, 3), p_bonferroni = signif(p_bonferroni, 3)), n = 50)

format_covariate_label <- function(term_names) {
label <- as.character(term_names)
label <- gsub("^Phases", "Phase ", label)
label <- gsub("^disease_category", "Disease category: ", label)
label <- gsub("^intervention_type", "Intervention: ", label)
label <- gsub("^log_enrollment$", "Log enrolment", label)
label <- gsub("^n_sites$", "Number of sites (per 10)", label)
label <- gsub("^is_multiregional$", "Multiregional conduct", label)
label <- gsub("^start_year$", "Start year", label)
label <- gsub("_", " ", label)
label <- gsub("PHASE ?", "Phase ", label)
label <- gsub("EARLY", "Early", label)
label <- gsub("\\|", "/", label)
trimws(label)
}

# Trial Characteristics
characteristics <- rq1_model_sponsorclusters %>%
filter(!grepl("^has_", term), term != "(Intercept)") %>%
mutate(scale_factor = ifelse(term == "n_sites", 10, 1),
percent_change = (exp(estimate * scale_factor) - 1) * 100,
percent_change_low = (exp(conf.low * scale_factor) - 1) * 100,
percent_change_high = (exp(conf.high * scale_factor) - 1) * 100)

print(characteristics %>% transmute(covariate = format_covariate_label(term),
percent_change = round(percent_change, 2),
percent_change_low = round(percent_change_low, 2),
percent_change_high = round(percent_change_high, 2),
p_value = signif(p.value, 3)), n = 40)

# RQ1 Association with T Duration
main_countries <- countries %>% filter(is_country)
country_order <- main_countries$country[order(main_countries$months)]
figure2_data <- main_countries %>% mutate(country = factor(country, levels = country_order))

figure2 <- ggplot(figure2_data, aes(months, country, colour = tier)) +
theme_minimal(base_size = 11) +
geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0, linewidth = 0.8, alpha = 0.65) +
geom_point(size = 2.6) +
  
scale_colour_manual(values = setNames(c("#2C7BB6", "#9E9E9E"), c("Significant after Bonferroni correction", "Not significant")), drop = FALSE, name = NULL) +
labs(x = "fewer months  ·  Months difference (95% CI)  ·  more months", y = NULL,
title = "Country Associations with Trial Duration",
subtitle = paste0("US baseline = 20.2 months |  Adjusted R² = ", round(rq1_model_summary$adj.r.squared, 3))) +
theme(legend.position = "bottom", panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
plot.title = element_text(face = "bold"),
plot.subtitle = element_text(colour = "grey35", size = 9))
ggsave(file.path(output_dir, "Figure2 - Country Associations.png"), figure2, width = 8.5, height = 6.5, dpi = 300)

# Table output
print(figure2_data %>% select(country, months, ci_low, ci_high, tier) %>% 
arrange(months) %>% mutate(months = round(months, 2), ci_low = round(ci_low, 2), ci_high = round(ci_high, 2)))

# Estimate vs Trial Volume
figure3 <- ggplot(main_countries, aes(n_trials, months, colour = tier)) +
theme_minimal(base_size = 11) +
geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, alpha = 0.45) +
geom_point(size = 2.5) +
geom_text_repel(aes(label = country), size = 2.8,
show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
seed = 42, box.padding = 0.45, point.padding = 0.25) +
scale_x_log10(breaks = c(1000, 1500, 2000, 3000, 4000, 5000),
minor_breaks = NULL, labels = scales::comma, expand = expansion(mult = c(0.08, 0.18))) +
scale_y_continuous(breaks = scales::breaks_pretty(7)) +
scale_colour_manual(values = setNames(c("#2C7BB6", "#9E9E9E"), c("Significant after Bonferroni correction", "Not significant")), drop = FALSE, name = NULL) +
labs(x = "Number of trials (log scale)", y = "fewer months  ·  Months difference (95% CI)  ·  more months",
title = "Country Estimate Precision vs Trial Count", subtitle = "US baseline = 20.2 months")+
theme(legend.position = "bottom",
panel.grid.minor.x = element_blank(),
panel.grid.major.y = element_line(colour = "grey92"),
panel.grid.major.x = element_line(colour = "grey95"),
plot.title = element_text(face = "bold"),
plot.subtitle = element_text(colour = "grey35", size = 9))
ggsave(file.path(output_dir, "Figure3 - Estimate Precision by Trial Volume.png"), figure3, width = 8.5, height = 6, dpi = 300)

# cook distance
cooks_d <- cooks.distance(rq1_model)
cooks_threshold <- 4 / nobs(rq1_model)
influential <- cooks_d > cooks_threshold

cat("Threshold =", round(cooks_threshold, 6), " Influential observations =", sum(influential),"% =", round(mean(influential) * 100, 2), "%")

png(file.path(output_dir, "Figure - RQ1 Cook's Distance.png"),
width = 3000, height = 2000, res = 300)

plot(cooks_d, type = "h",col = "black",
lwd = 1, xlab = "Observation index", ylab = "Cook's distance",
main = "RQ1 Model: Cook's Distance")

abline(h = cooks_threshold, col = "#D7301F", lty = 2,lwd = 1.5)
legend("topright",legend = paste0("Threshold = 4/n = ", format(cooks_threshold, digits = 6, scientific = FALSE)),
col = "#D7301F", lty = 2, lwd = 1.5, bty = "n", cex = 1)
dev.off()

trial_data_no_influential <- trial_data[!influential, ]; rq1_model_no_influential <- lm(rq1_formula, data = trial_data_no_influential)

print(tibble(model = c("Orig RQ1", "RQ1 without influential observations"),
n_trials = c(nobs(rq1_model), nobs(rq1_model_no_influential)),
R_squared = round(c(summary(rq1_model)$r.squared,summary(rq1_model_no_influential)$r.squared), 4),
adjusted_R_squared = round(c(summary(rq1_model)$adj.r.squared,summary(rq1_model_no_influential)$adj.r.squared), 4)))

# Calc of country shifts and significance after removing influential trials
countries_no_influential <- cluster_robust_estimates(rq1_model_no_influential, trial_data_no_influential) %>%
filter(term %in% rq1_countries) %>%
transmute(country = gsub("_", " ", gsub("has_", "", term)), no_influential_months = to_months(estimate, us_baseline),
no_influential_p = pmin(p.value * n_country_tests, 1))
cooks_comparison <- main_countries %>% transmute(country, primary_months = months, primary_p = p_bonferroni) %>%
left_join(countries_no_influential, by = "country") %>% mutate(shift_months = abs(no_influential_months - primary_months))
cat("Largest country shift after removing influential trials:", round(max(cooks_comparison$shift_months), 2), "months")
print(cooks_comparison %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = 50)

# Hold out validation
set.seed(42)
train_index  <- sample(seq_len(nrow(trial_data)), floor(0.8 * nrow(trial_data)))
train <- trial_data[train_index, ]
test  <- trial_data[-train_index, ] %>% filter(Phases %in% unique(train$Phases),
disease_category %in% unique(train$disease_category),
intervention_type %in% unique(train$intervention_type))
train_model  <- lm(rq1_formula, data = train)
predictions  <- predict(train_model, newdata = test)
r_squared_out_of_sample <- 1 - sum((test$log_duration - predictions)^2) / sum((test$log_duration - mean(test$log_duration))^2)

print(tibble(output = c("In sample R²", "Out of sample R²", "Median absolute error (months)"),
value = c(round(rq1_model_summary$r.squared, 4), round(r_squared_out_of_sample, 4),
round(median(abs(test$duration_months - exp(predictions))), 2))))

# Comparison, predicting every test trial at the training median
cat("Median absolute error at the training median:", round(median(abs(test$duration_months - median(train$duration_months))), 2), "months")

country_terms <- grep("^has_", names(coef(rq1_model)), value = TRUE)
coefficient_stability <- abs(to_months(coef(train_model)[country_terms], us_baseline) -
            to_months(coef(rq1_model)[country_terms], us_baseline))
print(tibble(largest_difference_months = round(max(coefficient_stability), 2),
countries_below_0.2_months = sum(coefficient_stability < 0.2),
total_countries = length(coefficient_stability)))

raw_scale_model <- lm(update(rq1_formula, duration_months ~ .), data = trial_data)
draw_qq <- function(fitted_model, panel_title) {
standardised_residuals <- rstandard(fitted_model)
qq_points <- qqnorm(standardised_residuals, plot.it = FALSE)
plot(qq_points$x, qq_points$y, pch = 16, cex = 0.2, col = adjustcolor("grey30", alpha.f = 0.25),
xlab = "Quantiles", ylab = "Standardised residuals", main = panel_title)
qqline(standardised_residuals, col = "#D7301F", lwd = 1.5)
}

# 4 plots, residuals vs fitted and normal Q-Q
draw_diagnostics <- function() {
par(mfrow = c(2, 2), mar = c(4.4, 4.4, 3, 1))
plot(fitted(raw_scale_model), resid(raw_scale_model), pch = 16, cex = 0.2,
col = adjustcolor("grey30", alpha.f = 0.25),
xlab = "Fitted values", ylab = "Residuals",
main = "(Raw duration) Residuals vs Fitted")
abline(h = 0, col = "#D7301F", lty = 2, lwd = 1.5)
draw_qq(raw_scale_model, "(Raw duration) Normal Q-Q")

plot(fitted(rq1_model), resid(rq1_model), pch = 16, cex = 0.2,
col = adjustcolor("grey30", alpha.f = 0.25),
xlab = "Fitted values", ylab = "Residuals",
main = "(Log duration) Residuals vs Fitted")
abline(h = 0, col = "#D7301F", lty = 2, lwd = 1.5)
draw_qq(rq1_model, "(Log duration) Normal Q-Q")
}
png(file.path(output_dir, "Figure5 - Model Diagnostics.png"), width = 3000, height = 2600, res = 300)
draw_diagnostics()
dev.off()

# RQ2 (multi-country trials)
multi_country_data <- trial_data %>% filter(n_countries >= 2)
print(count(multi_country_data))

model_A <- lm(log_duration ~ Phases + log_enrollment + disease_category + intervention_type + n_sites + start_year + n_countries + n_regions, data = multi_country_data) # no country added on
model_B <- lm(as.formula(paste("log_duration ~ Phases + log_enrollment + disease_category +",
    "intervention_type + n_sites + start_year +", "n_countries + n_regions +", country_string, "+ has_Other_country")), data = multi_country_data) # adds country_string - the indicators

adj_r_squared_model_a <- summary(model_A)$adj.r.squared; adj_r_squared_model_b <- summary(model_B)$adj.r.squared
partial_f_test <- anova(model_A, model_B)

print(tibble(model = c("RQ2 Model A (counts only)", "RQ2 Model B (counts + country identity)"),
n_trials = nrow(multi_country_data),
adj_r_squared = round(c(adj_r_squared_model_a, adj_r_squared_model_b), 5),
adj_r_squared_increase = c(NA, round(adj_r_squared_model_b - adj_r_squared_model_a, 5)),
partial_F = c(NA, round(partial_f_test$F[2], 2)),
p_value = c(NA, signif(partial_f_test$`Pr(>F)`[2], 3))))

# Multicollinearity Check
model_vif <- function(fitted_model) {vif_values <- vif(fitted_model)
if (is.matrix(vif_values)) vif_values[, "GVIF^(1/(2*Df))"]^2 else vif_values} #GVIF supports +2 levels of categories reason why this is used than vif plain formula
print(model_vif(rq1_model)); print(model_vif(model_A));print(model_vif(model_B))

# Country estimates in Model B when n_countries removed (VIF)
model_B_without_n_countries <- lm(update(formula(model_B), . ~ . - n_countries), data = multi_country_data)
vif_country_shift <- abs(to_months(coef(model_B_without_n_countries)[rq1_countries], us_baseline) - to_months(coef(model_B)[rq1_countries], us_baseline))
cat("Largest country shift with n_countries removed from Model B:", round(max(vif_country_shift), 2), "months")

# Sensitivity analyses

# RQ1 model as a function 
fit_countries <- function(data, country_list, label, baseline = us_baseline) {
pooled_columns <- setdiff(all_country_columns, c(country_list, "has_United_States"))
data$has_Other_country <- as.integer(rowSums(data[, pooled_columns, drop = FALSE]) > 0)
fitted_model <- lm(as.formula(paste("log_duration ~ Phases + log_enrollment + disease_category + intervention_type +",
"n_sites + is_multiregional + start_year +", paste(country_list, collapse = " + "),
"+ has_Other_country")), data = data)
cluster_robust_estimates(fitted_model, data) %>%
filter(grepl("^has_", term)) %>%
mutate(country = gsub("_", " ", gsub("has_", "", term)), analysis = label,
months = to_months(estimate, baseline),
ci_low = to_months(conf.low, baseline), ci_high = to_months(conf.high, baseline),
p_bonferroni = pmin(p.value * n_country_tests, 1),
percent = to_percent(estimate),
percent_low = to_percent(conf.low), percent_high = to_percent(conf.high),
adj_r_squared = summary(fitted_model)$adj.r.squared)
}

# Sensitivity analyses

# Letting th US as its own indicator
us_indicator_sensitivity <- fit_countries(trial_data, included_country_columns, "US indicator included")

cat("US indicator before R2", round(rq1_model_summary$adj.r.squared, 4), "After R2",
    round(us_indicator_sensitivity$adj_r_squared[1], 4))
print(us_indicator_sensitivity %>%
transmute(country, months = round(months, 2), ci_low = round(ci_low, 2),
ci_high = round(ci_high, 2), p_bonferroni = signif(p_bonferroni, 3)), n = Inf)

# US involvement across the data
us_coparticipation <- tibble(term = rq1_countries) %>%
mutate(country = gsub("_", " ", gsub("has_", "", term)),
percentage_with_US = sapply(term, function(each_term) 100 * mean(trial_data$has_United_States[trial_data[[each_term]] == 1])))

cat("US participation across all trials", round(100 * mean(trial_data$has_United_States), 1),"%")
print(us_coparticipation %>% select(country, percentage_with_US) %>% mutate(percentage_with_US = round(percentage_with_US, 1)),n = 150)

# germany omitted instead of the US, scaled at the same US baseline (20.2 months) so both specifications are comparable
germany_baseline_sensitivity <- fit_countries(trial_data, setdiff(included_country_columns, "has_Germany"), "Germany baseline", us_baseline)

# Temporal split
pre_2020_sensitivity  <- fit_countries(filter(trial_data, start_year <= 2019), rq1_countries, "Pre-COVID", us_baseline)
from_2020_sensitivity <- fit_countries(filter(trial_data, start_year >= 2020), rq1_countries, "COVID period", us_baseline)

# Pre vs Post Covid difference
period_shift <- pre_2020_sensitivity %>% filter(country != "Other country") %>%
transmute(country, before_2020_months = months,before_2020_percent = percent) %>%
left_join(from_2020_sensitivity %>% transmute(country, from_2020_months = months,from_2020_percent = percent),by = "country") %>%
mutate(shift_percentage_points = from_2020_percent - before_2020_percent) %>%
arrange(desc(abs(before_2020_months - from_2020_months))) %>%
  mutate(across(c(before_2020_months, from_2020_months), ~ round(.x, 2)), across(c(before_2020_percent, from_2020_percent, shift_percentage_points),~ round(.x, 1)))
print(period_shift, n = Inf)

# Lower threshold to 30
lower_threshold_cols <- setdiff(names(trial_counts[trial_counts >= 30]), "has_United_States")
lower_threshold_tidy <- fit_countries(trial_data, lower_threshold_cols, "30 trial threshold")
threshold_table <- countries %>% filter(is_country) %>% transmute(country, primary_est = months) %>%
  left_join(lower_threshold_tidy %>% transmute(country, low_thresh = months), by = "country") %>%
  mutate(diff = low_thresh - primary_est)
countries_count_under_trialthreshold <- setdiff(names(trial_counts[trial_counts >= 30]), "has_United_States"); length(countries_count_under_trialthreshold)

cat("Threshold of 30 -", length(lower_threshold_cols), "countries accounted for ",
sum(abs(threshold_table$diff) < 0.5, na.rm = TRUE),"of", nrow(threshold_table), "Move under 0.5 month ",
    "largest move", round(max(abs(threshold_table$diff), na.rm = TRUE), 2), "months")
print(threshold_table %>% transmute(country, primary_model_months = round(primary_est, 2),
threshold_30_months = round(low_thresh, 2), difference_months = round(diff, 2)),n = 100)

# Inclusion Treshold Comparison
threshold_comparison <- threshold_table %>% filter(!is.na(low_thresh)) %>%
pivot_longer(c(primary_est, low_thresh), names_to = "specification", values_to = "months") %>%
mutate(specification = factor(specification,levels = c("primary_est", "low_thresh"),
                                labels = c("(Within country model) 1,000 trial threshold", "30 trial threshold")),
country = factor(country, levels = country_order))

figure8 <- ggplot(threshold_comparison, aes(months, country, colour = specification, shape = specification)) +
theme_minimal(base_size = 11) +
geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
geom_line(aes(group = country), colour = "grey80", linewidth = 0.4, orientation = "y") +
geom_point(size = 2.3) +
scale_colour_manual(values = c("#2C7BB6", "#D7301F"), name = NULL) +
scale_shape_manual(values = c(16, 17), name = NULL) +
labs(x = "fewer months  ·  Months difference  ·  more months",
y = NULL, title = "Country effect vs inclusion threshold",
subtitle = "US baseline = 20.2 months") +  theme(legend.position = "bottom",
panel.grid.major.y = element_blank(),
plot.title = element_text(face = "bold"),
plot.subtitle = element_text(colour = "grey35", size = 9))
ggsave(file.path(output_dir, "Figure8 Inclusion Threshold.png"),figure8, width = 8.5, height = 6.5, dpi = 300)

# US vs Germany
countries_in_both_baselines <- intersect(main_countries$country, germany_baseline_sensitivity$country)
reference_comparison <- bind_rows(
  main_countries    %>% transmute(country, months, ci_low, ci_high, specification = "US omitted"),
  germany_baseline_sensitivity %>% transmute(country, months, ci_low, ci_high, specification = "Germany omitted")) %>%
  filter(country %in% countries_in_both_baselines) %>%
  mutate(country = factor(country, levels = country_order),
  specification = factor(specification, levels = c("US omitted", "Germany omitted")))

figure4_reference <- ggplot(reference_comparison, aes(months, country, colour = specification, shape = specification)) +
theme_minimal(base_size = 11) +
geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
position = position_dodge(width = 0.65), width = 0.22, linewidth = 0.5) +

geom_point(position = position_dodge(width = 0.65), size = 2.1) +
scale_colour_manual(values = c("US omitted" = "#2C7BB6", "Germany omitted" = "#D7301F"), name = NULL) +
scale_shape_manual(values = c(16, 15), name = NULL) +

labs(x = "fewer months  ·  Months difference (95% CI)  ·  more months", y = NULL,
       title = "Country Estimates with the US or Germany Indicator Omitted",
       subtitle = "US baseline = 20.2 months (both specifications)") +
theme(legend.position = "bottom",
panel.grid.major.y = element_blank(),
panel.grid.minor.x = element_blank(),
plot.title = element_text(face = "bold"),
plot.subtitle = element_text(colour = "grey35", size = 9))

ggsave(file.path(output_dir, "Baseline Specification Comparison.png"), figure4_reference, width = 8.5, height = 6.5, dpi = 300)

#table outputs
print(tibble(specification = c("Primary model (US omitted)", "Germany omitted"),
adjusted_R2 = round(c(summary(rq1_model)$adj.r.squared,
germany_baseline_sensitivity$adj_r_squared[1]), 4),
significant_countries = c(sum(main_countries$p_bonferroni < 0.05), sum(germany_baseline_sensitivity$p_bonferroni < 0.05)),
US_estimate_months = c(NA, round(germany_baseline_sensitivity$months[germany_baseline_sensitivity$country == "United States"], 1))))

# Figure 6 
period_comparison <- bind_rows(pre_2020_sensitivity, from_2020_sensitivity) %>% filter(country != "Other country") %>%
mutate(country  = factor(country, levels = country_order), analysis = factor(analysis, levels = c("Pre-COVID", "COVID period"),
labels = c("Pre-COVID (before 2020)", "COVID period (2020 onwards)")))

figure6_period <- ggplot(period_comparison, aes(months, country, colour = analysis, shape = analysis)) +
theme_minimal(base_size = 11) +
geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
position = position_dodge(width = 0.65), width = 0.22, linewidth = 0.5) +
geom_point(position = position_dodge(width = 0.65), size = 2.1) +
scale_colour_manual(values = c("Pre-COVID (before 2020)"     = "#2C7BB6", "COVID period (2020 onwards)" = "#D7301F"), name = NULL) +
scale_shape_manual(values = c(16, 15), name = NULL) +
labs(x = "fewer months  ·  Months difference (95% CI)  ·  more months", y = NULL,
title = "Country Associations with Trial Duration: Pre-COVID vs COVID Period",
subtitle = paste0("US baseline = 20.2 months", "")) +

theme(legend.position = "bottom",
panel.grid.major.y = element_blank(),
panel.grid.minor.x = element_blank(),
plot.title = element_text(face = "bold"),
plot.subtitle = element_text(colour = "grey35", size = 9))
ggsave(file.path(output_dir, "Figure6 Country Associations by Study Period.png"), figure6_period, width = 8.5, height = 6.5, dpi = 300)

# Sponsor adjustment (country coefficients estimated considering within sponsor variation)
trial_data_fixed_effects <- trial_data[as.numeric(table(trial_data$Sponsor)[trial_data$Sponsor]) >= 2, ]

predictor_matrix <- model.matrix(rq1_formula, data = trial_data_fixed_effects)[, -1, drop = FALSE]
adjusted_data <- as.data.frame(predictor_matrix)
adjusted_data$log_duration <- trial_data_fixed_effects$log_duration
adjusted_data$Sponsor      <- trial_data_fixed_effects$Sponsor

# subtract each sponsor's own average from every column, so only the differences within a sponsor's own trials are left
adjusted_data <- adjusted_data %>% group_by(Sponsor) %>%  mutate(across(where(is.numeric), ~ . - mean(.))) %>% ungroup()

# columns with no variation left (constant within every sponsor) can't be estimated
varying_columns <- names(adjusted_data)[
sapply(adjusted_data, function(x) is.numeric(x) && sd(x) > 1e-10)]; varying_columns <- setdiff(varying_columns, "log_duration")
fixed_effects_model <- lm(log_duration ~ . - 1, data = adjusted_data[, c("log_duration", varying_columns)])

n_clusters <- length(unique(trial_data_fixed_effects$Sponsor))
n_obs      <- nobs(fixed_effects_model)
n_params   <- length(coef(fixed_effects_model))
raw_vcov   <- vcovCL(fixed_effects_model, cluster = trial_data_fixed_effects$Sponsor, type = "HC0")
degrees_of_freedom_correction <- (n_clusters / (n_clusters - 1)) * ((n_obs - 1) / (n_obs - n_params - n_clusters)) #http://www.jstor.org/stable/24735989 pg 16 is formula, n_clusters subtracted as Sponsor FEs use up df
fixed_effects_vcov <- raw_vcov * degrees_of_freedom_correction

fixed_effects_results <- coeftest(fixed_effects_model, vcov. = fixed_effects_vcov, df = n_clusters - 1)
fixed_effects_table <- as.data.frame(fixed_effects_results[, 1:4]) %>%
rownames_to_column("term") %>% rename(estimate = Estimate, p.value = `Pr(>|t|)`) %>%
filter(grepl("^has_", term)) %>%
mutate(country = gsub("_", " ", gsub("has_", "", term)),
months  = to_months(estimate, us_baseline),
p_bonferroni = pmin(p.value * n_country_tests, 1))

r_squared_fixed_effects     <- 1 - sum(residuals(fixed_effects_model)^2) /
  sum((trial_data_fixed_effects$log_duration - mean(trial_data_fixed_effects$log_duration))^2)
adj_r_squared_fixed_effects <- 1 - (1 - r_squared_fixed_effects) * (n_obs - 1) / (n_obs - n_params - n_clusters)

# Sponsor trial volume on RQ1
trial_data$sponsor_volume <- log(as.numeric(table(trial_data$Sponsor)[trial_data$Sponsor]))
sponsor_volume_model <- lm(update(rq1_formula, . ~ . + sponsor_volume), data = trial_data)
country_shift <- abs(to_months(coef(sponsor_volume_model)[rq1_countries], us_baseline) - to_months(coef(rq1_model)[rq1_countries], us_baseline))
print(tibble(analysis = "Sponsor trial volume",adj_R2_primary = round(rq1_model_summary$adj.r.squared, 4),
             adj_R2_with_sponsor_volume = round(summary(sponsor_volume_model)$adj.r.squared, 4), largest_country_shift_months = round(max(country_shift), 2)))

primary_on_fe_sample <- lm(rq1_formula, data = trial_data_fixed_effects)

model_fit_comparison <- tibble(
model = c("Primary (RQ1)", "Primary (RQ1) on fixed-effects sample", "Sponsor fixed effects"),
n_trials = c(nrow(trial_data), n_obs, n_obs),
r_squared = c(round(rq1_model_summary$r.squared, 4),
round(summary(primary_on_fe_sample)$r.squared, 4),
round(r_squared_fixed_effects, 4)),
adj_r_squared = c(round(rq1_model_summary$adj.r.squared, 4),
round(summary(primary_on_fe_sample)$adj.r.squared, 4),
round(adj_r_squared_fixed_effects, 4)))
print(model_fit_comparison)

# Sponsor fixed effect RQ1 Comparison table
sponsor_ladder <- countries %>% filter(is_country) %>%
transmute(country, primary_months = round(months, 2), primary_p = round(p_bonferroni, 4)) %>%
left_join(fixed_effects_table %>% transmute(country, sponsor_adjusted_months = round(months, 2),
sponsor_adjusted_p = round(p_bonferroni, 4)),by = "country") %>% arrange(primary_months)
print(sponsor_ladder, n = 100)

# Table 6 country model comparison of with without fixed effects
fe_t_critical <- qt(0.975, df = n_clusters - 1)
table6 <- cluster_robust_estimates(primary_on_fe_sample, trial_data_fixed_effects) %>%
filter(term %in% rq1_countries) %>%
transmute(country = gsub("_", " ", gsub("has_", "", term)), country_model_months = to_months(estimate, us_baseline),
country_model_ci_low = to_months(conf.low, us_baseline), country_model_ci_high = to_months(conf.high, us_baseline)) %>%
left_join(fixed_effects_table %>% transmute(country, sponsor_fe_months = months,
sponsor_fe_ci_low = to_months(estimate - fe_t_critical * `Std. Error`, us_baseline),
sponsor_fe_ci_high = to_months(estimate + fe_t_critical * `Std. Error`, us_baseline),
sponsor_fe_p_bonferroni = p_bonferroni), by = "country") %>% arrange(country_model_months)
print(table6 %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = 500,  width = Inf)

# RQ2 Models A and B under sponsor fixed effects
multi_country_fixed_effects <- multi_country_data[as.numeric(table(multi_country_data$Sponsor)[multi_country_data$Sponsor]) >= 2, ]
n_sponsors_rq2 <- length(unique(multi_country_fixed_effects$Sponsor))

# Model A FE
predictor_matrix_A <- model.matrix(~ Phases + log_enrollment + disease_category + intervention_type + n_sites + start_year + n_countries + n_regions, data = multi_country_fixed_effects)[, -1, drop = FALSE]

adjusted_A <- as.data.frame(predictor_matrix_A)
adjusted_A$log_duration <- multi_country_fixed_effects$log_duration
adjusted_A$Sponsor      <- multi_country_fixed_effects$Sponsor
adjusted_A <- adjusted_A %>% group_by(Sponsor) %>%
  mutate(across(where(is.numeric), ~ . - mean(.))) %>% ungroup()
varying_A  <- names(adjusted_A)[sapply(adjusted_A, function(x) is.numeric(x) && sd(x) > 1e-10)]
model_A_fixed_effects <- lm(log_duration ~ . - 1, data = adjusted_A[, varying_A])

r_squared_A <- 1 - sum(residuals(model_A_fixed_effects)^2) /
  sum((multi_country_fixed_effects$log_duration - mean(multi_country_fixed_effects$log_duration))^2)
adj_r_squared_A <- 1 - (1 - r_squared_A) * (nobs(model_A_fixed_effects) - 1) /
  (nobs(model_A_fixed_effects) - length(coef(model_A_fixed_effects)) - n_sponsors_rq2)

# Model B Sp Fixed Effects
predictor_matrix_B <- model.matrix(as.formula(paste("log_duration ~ Phases + log_enrollment + disease_category +","intervention_type + n_sites + start_year +", "n_countries + n_regions +",
    country_string, "+ has_Other_country")), data = multi_country_fixed_effects)[, -1, drop = FALSE]

adjusted_B <- as.data.frame(predictor_matrix_B)
adjusted_B$log_duration <- multi_country_fixed_effects$log_duration
adjusted_B$Sponsor      <- multi_country_fixed_effects$Sponsor
adjusted_B <- adjusted_B %>% group_by(Sponsor) %>%
  mutate(across(where(is.numeric), ~ . - mean(.))) %>% ungroup()
varying_B  <- names(adjusted_B)[sapply(adjusted_B, function(x) is.numeric(x) && sd(x) > 1e-10)]
model_B_fixed_effects <- lm(log_duration ~ . - 1, data = adjusted_B[, varying_B])

r_squared_B <- 1 - sum(residuals(model_B_fixed_effects)^2) /
  sum((multi_country_fixed_effects$log_duration - mean(multi_country_fixed_effects$log_duration))^2)
adj_r_squared_B <- 1 - (1 - r_squared_B) * (nobs(model_B_fixed_effects) - 1) /
  (nobs(model_B_fixed_effects) - length(coef(model_B_fixed_effects)) - n_sponsors_rq2)

print(tibble(output = "RQ2 sponsor fixed effects", n_trials = nrow(multi_country_fixed_effects),
n_sponsors = n_sponsors_rq2, model_A_adj_R2 = round(adj_r_squared_A, 4), model_B_adj_R2 = round(adj_r_squared_B, 4)))

# RQ2 Models A and B on the same trials without sponsor fixed effects
model_A_same_trials <- lm(formula(model_A), data = multi_country_fixed_effects)
model_B_same_trials <- lm(formula(model_B), data = multi_country_fixed_effects)

# change in adjusted R-squared between Models A and B
r2_change_same_trials <- summary(model_B_same_trials)$adj.r.squared - summary(model_A_same_trials)$adj.r.squared

# change from the Sponsor fe models
r2_change_fixed_effects <- adj_r_squared_B - adj_r_squared_A

cat("Adjusted R-squared change without fixed effects",round(r2_change_same_trials, 4),
"Adjusted R-squared change with Sponsor fixed effects",round(r2_change_fixed_effects, 4))
cat("Complete")
