library(phylolm)
library(officer)
library(flextable)
library(dplyr)



# Extract and clean summary
summary_phylo <- summary(model_25)


coef_df <- as.data.frame(summary_phylo$coefficients) %>%
  rownames_to_column("Term") %>%
  rename(
    Estimate = Estimate,
    `Std. Error` = StdErr,
    `t value` = t.value,
    `p value` = p.value
  ) %>%
  mutate(
    Estimate = round(Estimate, 2),
    `Std. Error` = round(`Std. Error`, 2),
    `t value` = round(`t value`, 2),
    `p value` = formatC(`p value`, format = "e", digits = 2)  # scientific notation
  )

# Create flextable
ft <- flextable(coef_df) %>%
  autofit() %>%
  theme_vanilla() %>%
  set_caption("Phylogenetic Linear Model Coefficients") %>%
  align(align = "center", part = "all")

# Export to Word
doc <- read_docx() %>%
  body_add_par("Phylogenetic Linear Model Summary", style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = "phylolm_summary_pc2_woody.docx")
