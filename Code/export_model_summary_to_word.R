library(dplyr)
library(tibble)
library(flextable)
library(officer)
library(purrr)  # for map_chr

# Get model summary
summary_phylo <- summary(model_25)

# Rename terms
term_labels <- c(
  "(Intercept)" = "Intercept",
  "provenancenon_native" = "Provenance (non-native)",
  "ave_tmean" = "Mean Temp.",
  "ave_precip" = "Mean Precip.",
  "I(ave_tmean^2)" = "Mean Temp.^2",
  "I(ave_precip^2)" = "Mean Precip.^2",
  "provenancenon_native:ave_tmean" = "Provenance:Mean Temp.",
  "provenancenon_native:ave_precip" = "Provenance:Mean Precip."
)

# Build cleaned coefficient table with significance stars
coef_df <- as.data.frame(summary_phylo$coefficients) %>%
  rownames_to_column("Term_raw") %>%
  mutate(
    Term = map_chr(Term_raw, ~ term_labels[.x]),  # Correct mapping
    Estimate = round(Estimate, 2),
    `Std. Error` = round(StdErr, 2),
    `t value` = round(t.value, 2),
    `p value` = case_when(
      p.value <= 0.001 ~ paste0(formatC(p.value, format = "e", digits = 2), "***"),
      p.value <= 0.01  ~ paste0(formatC(p.value, format = "e", digits = 2), "**"),
      p.value <= 0.05  ~ paste0(formatC(p.value, format = "e", digits = 2), "*"),
      TRUE             ~ formatC(p.value, format = "e", digits = 2)
    ),
    Significant = p.value <= 0.05
  )

# Create flextable with bolded significant terms
ft <- flextable(coef_df, col_keys = c("Term", "Estimate", "Std. Error", "t value", "p value"))

# Apply bold formatting to significant terms row-by-row
for (i in seq_len(nrow(coef_df))) {
  ft <- compose(
    ft,
    i = i,
    j = "Term",
    value = as_paragraph(
      as_chunk(coef_df$Term[i], props = fp_text(bold = coef_df$Significant[i]))
    )
  )
}

# Final formatting
ft <- ft %>%
  autofit() %>%
  theme_vanilla() %>%
  set_caption("Phylogenetic Linear Model Coefficients") %>%
  align(j = "Term", align = "left", part = "all") %>%
  align(j = c("Estimate", "Std. Error", "t value", "p value"), align = "center", part = "all")

# Export to Word
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "phylolm_summary_PC2_Herbs.docx")

