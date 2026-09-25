rm(list = ls()); library(tidyverse);library(lubridate)
data <- read.csv("C:/Users/L042952/Desktop/Complete_Data.csv", header = TRUE)

# Region lists [https://unstats.un.org/unsd/methodology/m49/] my categorisation is simplified to my company 4 regions
NA_Countries <- c("United States", "Canada")

LATAM_Countries <- c(
  "Mexico", "Brazil", "Argentina", "Chile", "Colombia", "Peru",
  "Ecuador", "Guatemala", "Costa Rica", "Panama",
  "Dominican Republic", "Honduras", "El Salvador", "Uruguay",
  "Paraguay", "Bolivia", "Venezuela", "Nicaragua",
  "Trinidad and Tobago", "Jamaica", "Cuba", "Haiti", "Barbados",
  "Bahamas", "Suriname", "Guyana", "Belize")

EMEA_Countries <- c(
  "United Kingdom", "Germany", "France", "Spain", "Italy",
  "Netherlands", "Belgium", "Switzerland", "Austria", "Sweden",
  "Denmark", "Norway", "Finland", "Ireland", "Portugal", "Poland",
  "Czech Republic", "Hungary", "Romania", "Bulgaria", "Greece",
  "Croatia", "Slovakia", "Slovenia", "Lithuania", "Latvia",
  "Estonia", "Serbia", "Ukraine", "Russia", "Turkey", "Israel",
  "Saudi Arabia", "South Africa", "Egypt", "United Arab Emirates",
  "Lebanon", "Jordan", "Kenya", "Nigeria", "Ghana", "Morocco",
  "Tunisia", "Iran", "Iraq", "Kuwait", "Qatar", "Bahrain", "Oman",
  "Cyprus", "Malta", "Iceland", "Luxembourg", "North Macedonia",
  "Bosnia and Herzegovina", "Montenegro", "Albania", "Moldova",
  "Georgia", "Armenia", "Azerbaijan", "Ethiopia", "Tanzania",
  "Uganda", "Mozambique", "Zambia", "Zimbabwe", "Senegal",
  "Cameroon", "Ivory Coast", "Belarus", "Kazakhstan", "Algeria",
  "Democratic Republic of the Congo", "Burkina Faso", "Mali",
  "Gabon", "Rwanda", "Malawi", "Botswana", "Sudan", "The Gambia",
  "Central African Republic", "Sierra Leone", "Guinea",
  "Equatorial Guinea", "South Sudan", "Libya", "Mauritius",
  "Uzbekistan", "Monaco")

ASIAPAC_Countries <- c(
  "China", "Japan", "South Korea", "India", "Australia",
  "New Zealand", "Taiwan", "Hong Kong", "Singapore", "Malaysia",
  "Thailand", "Philippines", "Indonesia", "Vietnam", "Pakistan",
  "Bangladesh", "Sri Lanka", "Myanmar", "Cambodia", "Nepal",
  "Mongolia")

all_known <- c(NA_Countries, LATAM_Countries,
               EMEA_Countries, ASIAPAC_Countries)

# Standardisation
country_aliases <- c(
  "Turkey (Türkiye)" = "Turkey",
  "Czechia" = "Czech Republic",
  "Côte d'Ivoire" = "Ivory Coast",
  "Korea, Republic of" = "South Korea",
  "Republic of Korea" = "South Korea",
  "Russian Federation" = "Russia",
  "Serbia and Montenegro" = "Serbia",
  "Martinique" = "France",
  "Guadeloupe" = "France",
  "USA" = "United States",
  "United States Minor Outlying Islands" = "United States",
  "Puerto Rico" = "United States")

standardise_country <- function(name) {
  if (length(name) == 0 || is.na(name) || name == "") return(NA_character_)
  if (name %in% names(country_aliases)) return(country_aliases[[name]])
  return(name)
}

# Disease categorisation from MeSH 2026
mesh_raw <- readLines("C:/Users/L042952/Desktop/mesh_terms.bin")
mesh_raw <- trimws(mesh_raw)
mesh_raw <- mesh_raw[mesh_raw != ""] # removes empty lines

# Build the mesh lookup table one line at a time
mesh_terms <- character(length(mesh_raw))
mesh_codes <- character(length(mesh_raw))
for (i in seq_along(mesh_raw)) {
  parts <- strsplit(mesh_raw[i], ";")[[1]]
  mesh_terms[i] <- trimws(parts[1])
  mesh_codes[i] <- trimws(parts[2])
}

mesh <- data.frame(term = tolower(mesh_terms), code = mesh_codes, stringsAsFactors = FALSE)
mesh$chapter <- substr(mesh$code, 1, 3) # extracts first 3 characters of each MeSH code

# Chapter to 16 categories - # visible in https://meshb.nlm.nih.gov/treeView under Diseases [C]
category_mapping <- list(
  Oncology         = c("C04"),
  Cardiometabolic  = c("C14", "C18", "C19"),
  Neuro_Psych      = c("C10", "F03"),
  Infectious       = c("C01"),
  Immune           = c("C20"),
  Respiratory      = c("C08"),
  Musculoskeletal  = c("C05"),
  Digestive        = c("C06"),
  Urogenital_Repro = c("C12"),
  Dermatology      = c("C17"),
  Hematology       = c("C15"),
  Ophthalmology    = c("C11"),
  ENT              = c("C09"),
  Dental           = c("C07"),
  Congenital_Rare  = c("C16"),
  Injury           = c("C26"))

# Creates the terms for each category, one at a time
patterns <- list()
for (category_name in names(category_mapping)) {
  chapters_wanted <- category_mapping[[category_name]]
  matching_terms <- mesh$term[mesh$chapter %in% chapters_wanted & nchar(mesh$term) >= 3]
  patterns[[category_name]] <- paste(matching_terms, collapse = "|")
}
patterns[["Oncology"]] <- paste(patterns[["Oncology"]], "cancer", "tumour", "tumor", sep = "|") # resolves for oncology freq terms

# Sets the order for assigning disease categories, from top to bottom, it is first match approach, while diseases have lot of cross-over in categories - given use of this information - this approach is fine
category_priority <- c(
  "Oncology", "Congenital_Rare", "Infectious", "Immune",
  "Neuro_Psych", "Cardiometabolic", "Respiratory", "Digestive",
  "Urogenital_Repro", "Hematology", "Dermatology", "Ophthalmology",
  "Musculoskeletal", "ENT", "Dental", "Injury")

# I identified location was like 'Boston, Massachusetts, United States' so I split it into parts to find the country
find_country_in_parts <- function(parts) {
  number_of_parts <- length(parts)
  
last_part <- parts[number_of_parts]
standardised_name <- standardise_country(last_part)
if (!is.na(standardised_name) && standardised_name %in% all_known) return(standardised_name)
  
if (number_of_parts >= 2) {
second_last_part <- parts[number_of_parts - 1]
standardised_name <- standardise_country(second_last_part)
if (!is.na(standardised_name) && standardised_name %in% all_known) return(standardised_name)
}
# parts means say we have Maradyke, Cork, Ireland (3) versus Cork, Ireland
if (number_of_parts >= 3) {
third_last_part <- parts[number_of_parts - 2]
standardised_name <- standardise_country(third_last_part)
if (!is.na(standardised_name) && standardised_name %in% all_known) return(standardised_name)
}
  
# no country found among the last three parts
  return(NA_character_)
}

# Extract countries from Locations column
get_countries <- function(location_string) {
  if (is.na(location_string) || location_string == "") return(character(0))
  sites <- unlist(strsplit(as.character(location_string), "\\|"))
  countries <- c()
  for (site in sites) {
    parts <- trimws(unlist(strsplit(site, ",")))
    parts <- parts[parts != ""]
    if (length(parts) == 0) next
    found_country <- find_country_in_parts(parts)
    if (!is.na(found_country)) countries <- c(countries, found_country)
  }
  return(unique(countries))
}

# Assign region
get_region <- function(country) {
  if (country %in% NA_Countries)      return("NA")
  if (country %in% LATAM_Countries)   return("LATAM")
  if (country %in% EMEA_Countries)    return("EMEA")
  if (country %in% ASIAPAC_Countries) return("ASIAPAC")
  return("Other")
}

# Parse dates, this is a clean up
data$start_dt <- ymd(data$Start.Date, quiet = TRUE)
data$primary_dt <- ymd(data$Primary.Completion.Date, quiet = TRUE)

missing_start_dates <- is.na(data$start_dt)
missing_primary_dates <- is.na(data$primary_dt)
data$start_dt[missing_start_dates] <- ym(data$Start.Date[missing_start_dates], quiet = TRUE)
data$primary_dt[missing_primary_dates] <- ym(data$Primary.Completion.Date[missing_primary_dates], quiet = TRUE)

# Calculate variables
data$duration_months <- as.numeric(difftime(data$primary_dt, data$start_dt, units = "days")) / 30.4375
data$start_year <- year(data$start_dt)

data$countries <- lapply(data$Locations, get_countries)
data$n_countries <- sapply(data$countries, length)

data$n_sites <- str_count(data$Locations, "\\|") + 1
data$n_sites[is.na(data$Locations) | data$Locations == ""] <- 0

data$regions <- lapply(data$countries, function(country_list) unique(sapply(country_list, get_region)))
data$n_regions <- sapply(data$regions, length)
data$has_NA <- sapply(data$regions, function(region_list) as.integer("NA" %in% region_list))
data$has_LATAM <- sapply(data$regions, function(region_list) as.integer("LATAM" %in% region_list))
data$has_EMEA <- sapply(data$regions, function(region_list) as.integer("EMEA" %in% region_list))
data$has_ASIAPAC <- sapply(data$regions, function(region_list) as.integer("ASIAPAC" %in% region_list))
data$is_multiregional <- as.integer(data$n_regions > 1)

# Converts text before the colon (taking advantage of a data formatting logic for example "DRUG: Paracetamol 500mg", isolated to "DRUG", the idea is to weight the model against different drug types
data$intervention_type <- word(data$Interventions, 1, sep = ":")

# Healthy volunteers involved in Ph1 trials - which are in scope - so accounted for in category
condition_lowercase <- tolower(data$Conditions)
disease_category <- rep("Other", length(condition_lowercase))
disease_category[grepl("healthy|normal volunteer", condition_lowercase)] <- "Healthy_Volunteer"
for (category_name in rev(category_priority)) {
  matches <- grepl(patterns[[category_name]], condition_lowercase)
  disease_category[matches] <- category_name
}
disease_category[is.na(condition_lowercase) | condition_lowercase == ""] <- "Unknown"
data$disease_category <- disease_category

data$Enrollment <- as.numeric(data$Enrollment)
data$log_enrollment <- log(data$Enrollment)

# Exclusion criteria
data <- data %>%
filter(grepl("^NCT", NCT.Number),
!is.na(start_dt),
!is.na(primary_dt),
Enrollment > 0,
duration_months > 0,
n_countries >= 1,
!grepl("PHASE4", Phases, ignore.case = TRUE))

# Binary country indicators
all_countries_list <- names(sort(table(unlist(data$countries)), decreasing = TRUE))
for (country in all_countries_list) {
  column_name <- paste0("has_", gsub(" ", "_", country))
  data[[column_name]] <- sapply(data$countries, function(country_list) as.integer(country %in% country_list))}

data$log_duration <- log(data$duration_months)

trial_data <- data %>% select(-countries, -regions, -Locations)
write.csv(trial_data, "C:/Users/L042952/Desktop/Clean_Data.csv", row.names = FALSE)
cat("Trials", nrow(trial_data), "Sponsors", length(unique(trial_data$Sponsor)),
    " Largest sponsor", max(table(trial_data$Sponsor)), "trials")
