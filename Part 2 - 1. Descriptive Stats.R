rm(list = ls()); library(tidyverse); library(car); library(ggrepel); library(sandwich); library(lmtest)
trial_data <- read.csv( "C:/Users/L042952/Desktop/Clean_Data.csv", header = TRUE, stringsAsFactors = FALSE)
output_dir  <- "C:/Users/L042952/Desktop"

# Descriptive stats table
total_trials <- nrow(trial_data)
print(table(trial_data$Phases));print(table(trial_data$intervention_type));print(table(trial_data$disease_category))

print(round(prop.table(table(trial_data$Phases)) * 100, 2))
print(round(prop.table(table(trial_data$intervention_type)) * 100, 2))
print(round(prop.table(table(trial_data$disease_category)) * 100, 2))

# Quantile function for calculation
get_quantile <- function(x, quantile_position) round(quantile(x, quantile_position), 1)

# Trial duration
cat( "Median Duration months", round(median(trial_data$duration_months), 1),
     "IQR", get_quantile(trial_data$duration_months, .25), "-", get_quantile(trial_data$duration_months, .75))

cat("Median Enrolment:", round(median(trial_data$Enrollment), 1),
    "IQR", get_quantile(trial_data$Enrollment, .25), "-", get_quantile(trial_data$Enrollment, .75), "\n")

cat("Median Sites:", round(median(trial_data$n_sites), 1), 
    "IQR", get_quantile(trial_data$n_sites, .25), "-", get_quantile(trial_data$n_sites, .75), "\n")

single_country <- sum(trial_data$n_countries == 1);  multi_country <- sum(trial_data$n_countries >= 2)
cat("Single country:", single_country,"(", round(single_country / total_trials * 100, 2), "%)",
    "Multi-country:", multi_country, "(", round(multi_country / total_trials * 100, 2), "%)")

pre_covid <- sum(trial_data$start_year <= 2019)
covid_post_covid <- sum(trial_data$start_year >= 2020 & trial_data$start_year <= 2025)

cat("Before 2020:", pre_covid, "(", round(pre_covid / total_trials * 100, 2), "%)",  " After 2020 until 2025:", covid_post_covid,
    "(", round(covid_post_covid / total_trials * 100, 2), "%)", "Earliest start year:", min(trial_data$start_year), "\n")