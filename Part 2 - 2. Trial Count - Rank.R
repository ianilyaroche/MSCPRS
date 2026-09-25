rm(list = ls()); library(tidyverse); library(car); library(ggrepel); library(sandwich); library(lmtest)
trial_data <- read.csv( "C:/Users/L042952/Desktop/Clean_Data.csv", header = TRUE, stringsAsFactors = FALSE)
output_dir  <- "C:/Users/L042952/Desktop"

# Removes the region columns from graph, its filtering by 'has_'
region_columns <- c("has_NA", "has_LATAM", "has_EMEA", "has_ASIAPAC")
all_country_columns <- setdiff(grep("^has_", names(trial_data), value = TRUE), region_columns)
trial_counts <- sort(sapply(all_country_columns, function(country_column) sum(trial_data[[country_column]])), decreasing = TRUE)

included_country_columns <- names(trial_counts[trial_counts >= 1000]) ; length(included_country_columns)
below_threshold_columns <- setdiff(all_country_columns, included_country_columns); length(below_threshold_columns)
trial_data$has_Other_country <- as.integer(rowSums(trial_data[, below_threshold_columns, drop = FALSE]) > 0)
print(trial_counts[trial_counts >= 1000])

# Trial count graph
draw_figure1 <- function() {
par(mar = c(4.4, 4.8, 3, 1.2)) 
ranked_counts  <- as.numeric(trial_counts)
entered_individually <- ranked_counts >= 1000
point_colours    <- ifelse(entered_individually, "#2C7BB6", "#9E9E9E")

plot(seq_along(ranked_counts), ranked_counts,
type = "n",
xlab = "Country rank",
ylab = "Number of trials",
main = "Number of trials vs Country rank")

lines(seq_along(ranked_counts), ranked_counts, col = "grey70", lwd = 0.7)
abline(h = 1000, lty = 2, col = "#D7301F")
points(seq_along(ranked_counts), ranked_counts, pch = 21, cex = 0.8,
col = point_colours, bg = point_colours)
text(x = 1, y = 1000, adj = c(0, -0.6), cex = 0.75, col = "grey35",
       
labels = paste0("                          Inclusion threshold = 1,000 trials"))
legend("topright", bty = "n", cex = 0.8, pch = 21,
col = c("#2C7BB6", "#9E9E9E"), pt.bg = c("#2C7BB6", "#9E9E9E"),
         legend = c(paste0("Entered individually (", length(included_country_columns), ")"),
                    "Pooled as 'Other country'"))
}
draw_figure1()
png(file.path(output_dir, "Figure1 Trial Counts.png"),width = 2400, height = 1500, res = 300)
draw_figure1()
dev.off() # ends and saves