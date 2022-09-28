library(ggplot2)
library(dplyr)


taxonomy_table <- read.csv("Cqui_family_100_all.EEs.tax.tsv", sep="\t", header = TRUE)
molecule_type <- taxonomy_table %>% count(Molecule_type)
taxonomy_table <- taxonomy_table[taxonomy_table$Molecule_type != "dsDNA-RT", ]
taxonomy_table <- taxonomy_table[taxonomy_table$Molecule_type != "ssRNA-RT", ]

family_count <- taxonomy_table %>% count(Family)
genus_count <- taxonomy_table %>% count(Genus)


##MOLECULE TYPE

# Compute percentages
molecule_type$fraction = molecule_type$n / sum(molecule_type$n)

# Compute the cumulative percentages (top of each rectangle)
molecule_type$ymax = cumsum(molecule_type$fraction)

# Compute the bottom of each rectangle
molecule_type$ymin = c(0, head(molecule_type$ymax, n=-1))

# Compute label position
molecule_type$labelPosition <- (molecule_type$ymax + molecule_type$ymin) / 2

# Compute a good label
molecule_type$label_value <- paste0("N = ", molecule_type$n)
# Make the plot
ggplot(molecule_type, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Molecule_type)) +
  geom_rect() +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
  geom_label( x=3.5, aes(y=labelPosition, label=label_value), size=6) +
  theme_void()

##FAMILY

# Compute percentages
family_count$fraction = family_count$n / sum(family_count$n)

# Compute the cumulative percentages (top of each rectangle)
family_count$ymax = cumsum(family_count$fraction)

# Compute the bottom of each rectangle
family_count$ymin = c(0, head(family_count$ymax, n=-1))

# Compute label position
family_count$labelPosition <- (family_count$ymax + family_count$ymin) / 2

# Compute a good label
family_count$label_value <- paste0("N = ", family_count$n)
# Make the plot
ggplot(family_count, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Family)) +
  geom_rect() +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
  geom_label( x=3.5, aes(y=labelPosition, label=label_value), size=6) +
  theme_void()
