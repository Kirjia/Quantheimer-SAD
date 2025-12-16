library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(janitor)
library(corrplot)
library(fastDummies)
library(hrbrthemes)
library(heatmaply)


dt <- read_tsv(file = 'C:/Users/enzo2/OneDrive/Desktop/openneuro/participants.tsv', na = c("", "NA", "n/a", "N/A"))

df<-data.frame(dt)
varianza <- data.frame()

df <- df %>%
  mutate(across(9:87 & where(is.character), readr::parse_number))


ggplot(data = df, aes(x = APOE_haplotype, fill = factor(APOE_haplotype))) + geom_bar() + theme_minimal() +
  labs(
    title = "Genes presence",
    ylab = "Count",
    fill = "Genes"
  )

ggplot(data = df, aes(x = APOE_rs429358, fill = factor(APOE_rs429358))) + geom_bar() + theme_minimal() +
  labs(
    title = "apoe rs429358 presence",
    ylab = "Count",
    fill = "Genes"
  )

ggplot(data = df, aes(x = APOE_rs7412, fill = factor(APOE_rs7412))) + geom_bar() + theme_minimal() +
  labs(
    title = "apoe rs7412",
    ylab = "Count",
    fill = "Genes"
  )

ggplot(data = df, aes(y = BDI)) + geom_boxplot() + labs(
  title = "BDI boxplot"
)

#density plot
df %>% ggplot( aes(x=BDI)) +
  geom_density(fill="#696089", color="#000000", alpha=0.8) +
  ggtitle("Distribuzione del Beck depression inventory") + theme_ipsum()

df %>% ggplot( aes(x=RPM)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione di RPM(Raven's progressive Matrix)") + theme_ipsum()

df %>% ggplot( aes(x=SES)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. rapporto socio economico SES") + theme_ipsum()

df %>% ggplot( aes(x=EHI)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. EHI grado di quanto una persona sia destrorsa") + theme_ipsum()

df %>% ggplot( aes(x=AUDIT)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. indice alcoolemico") + theme_ipsum()

df %>% ggplot( aes(x=education)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. educazione") + theme_ipsum()

df %>% ggplot( aes(x=BMI)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. massa corporea") + theme_ipsum()

df %>% ggplot( aes(x=learning_deficits)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. indice defict apprendimento") + theme_ipsum()



df %>% ggplot( aes(x=BDI, y=RPM)) +
  geom_point(size=6, color="#69b3a2") +
  ggtitle("Transparency")

df_long <- pivot_longer(df,
                        cols = c(NEO_NEU, NEO_EXT, NEO_OPE, NEO_AGR, NEO_CON),
                        names_to = "Variabile",
                        values_to = "Valore")

# 2. Ora plottiamo tutto insieme
ggplot(df_long, aes(x = Valore, fill = Variabile)) +
  geom_density(alpha = 0.6, color = NA) + # alpha rende trasparente
  facet_wrap(~ Variabile, ncol = 2) +     # <--- Crea la griglia 2x2
  theme_ipsum() +
  labs(title = "Confronto Densità dei valori NEO", fill = "Variabile")

df <- df %>%
  mutate(
    APOE_risk_score = case_match(APOE_haplotype,
                                 c("e2/e2", "e3/e2") ~ 1,  # Protettivo
                                 "e3/e3"             ~ 2,  # Neutro
                                 c("e3/e4", "e2/e4") ~ 3,  # Rischio 1
                                 "e4/e4"             ~ 4,  # Rischio 2
                                 .default = NA             # Gestione errori
    )
  )



df$PICALM_rs3851179 <- factor(df$PICALM_rs3851179)
df <- dummy_cols(df, select_columns = ("PICALM_rs3851179"), remove_first_dummy = TRUE )



boxplot(na.omit(df$CVLT_12), xlab="CVLT total hits", col = "green")
df %>% ggplot(aes(y = CVLT_12))+geom_boxplot( outlier.color = "black", outlier.shape = 16, outlier.size = 2, fill="red", ylab= "count")
boxplot(df$education)

df %>% ggplot(aes(y = df$RPM))+geom_boxplot( outlier.color = "black", outlier.shape = 16, outlier.size = 2, fill="red")



tmp <- df[9:90]
df_clean <- tmp[, apply(tmp, 2, var, na.rm = TRUE) != 0]
df_clean[is.na(df_clean)] <- 0


cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x),
    colnames(x),
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
p <- cor.test.p(df_clean)

heatmaply_cor(
  cor(df_clean),            # Questo dice a R di scriverli
  # Opzionale: Estetica del testo
  node_type = "scatter",
  hclust_method= "ward.D2",
  point_size_mat = -log10(p),
  point_size_name = "-log10(p-value)",
  dendrogram = "none",
  label_names = c("x", "y", "Correlation")
)

