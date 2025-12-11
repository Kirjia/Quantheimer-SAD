library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(janitor)
library(corrplot)

dt <- read.table(file = 'C:/Users/enzo2/OneDrive/Desktop/openneuro/participants.tsv', sep = '\t', header = TRUE)

df<-data.frame(dt)

ggplot(data = df, aes(x = APOE_haplotype, fill = factor(APOE_haplotype))) + geom_bar() + theme_minimal() +
  labs(
    title = "Genes presence",
    ylab = "Count",
    fill = "Genes"
  )

df <- df %>%
  mutate(across(
    9:87 & where(is.character),
    ~ parse_number(., locale = locale(decimal_mark = ",", grouping_mark = "."))
  ))

boxplot(na.omit(df$CVLT_12), xlab="CVLT total hits", col = "green")
df %>% ggplot(aes(y = CVLT_12))+geom_boxplot( outlier.color = "black", outlier.shape = 16, outlier.size = 2, fill="red", ylab= "count")
boxplot(df$education)

df %>% ggplot(aes(y = df$RPM))+geom_boxplot( outlier.color = "black", outlier.shape = 16, outlier.size = 2, fill="red")

df <- df %>%
  mutate(across(
    1:9 & where(is.character),
    ~ parse_number(., locale = locale(decimal_mark = ",", grouping_mark = "."))
  ))

corrplot(df)

