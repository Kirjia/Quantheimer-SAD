if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, ggthemes, dplyr, tidyr, readr, janitor, corrplot, fastDummies, heatmaply, factoextra, skimr, psych)




dt <- read_tsv(file = file.path(getwd(), "/dataset/participants.tsv"), na = c("", "NA", "n/a", "N/A"))

df<-data.frame(dt)

head(df)

#quanti NA ci sono in tutto il dataset
sum(is.na(df))

#quali colonne hanno NA
colSums(is.na(df))

df_stats <- skim(df)
df_py_stats <- describe(df)

df <- df %>%
  mutate(across(9:87 & where(is.character), readr::parse_number))
df$APOE_rs429358[122] = "T/T"

df <- df %>% slice(-108)


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

ggplot(data = df, aes(x= diabetes, fill = factor(diabetes))) + geom_bar() + theme_minimal() +
  labs(
    title = "distribuzione del diabete",
    ylab = "Count",
    fill= "diabetes"
  )

ggplot(data = df, aes(x= hypertension, fill = factor(hypertension))) + geom_bar() + theme_minimal() +
  labs(
    title = "distribuzione hypertension",
    ylab = "Count",
    fill= "hypertension"
  )

ggplot(data = df, aes(y = BDI)) + geom_boxplot() + labs(
  title = "BDI boxplot"
)

#density plot
df %>% ggplot( aes(x=BDI)) +
  geom_density(fill="#696089", color="#000000", alpha=0.8) +
  ggtitle("Distribuzione del Beck depression inventory") + theme_minimal()

df %>% ggplot( aes(x=RPM)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione di RPM(Raven's progressive Matrix)") + theme_minimal()

df %>% ggplot( aes(x=SES)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. rapporto socio economico SES") + theme_minimal()

df %>% ggplot( aes(x=EHI)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. EHI grado di quanto una persona sia destrorsa") + theme_minimal()

df %>% ggplot( aes(x=AUDIT)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. indice alcoolemico") + theme_minimal()

df %>% ggplot( aes(x=education)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. educazione") + theme_minimal()

df %>% ggplot( aes(x=BMI)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. massa corporea") + theme_minimal()

df %>% ggplot( aes(x=learning_deficits)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distr. indice defict apprendimento") + theme_minimal()



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
  theme_minimal() +
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

df<- df%>%
  mutate(
    APOE_rs429358 = case_match(APOE_rs429358,
                               c("T/T") ~ 1,
                               c("T/C") ~ 2,
                               c("C/T") ~ 3,
                               c("C/C") ~ 4,
                               .default = NA)
  )

df<- df%>%
  mutate(
    APOE_rs7412 = case_match(APOE_rs7412,
                               c("T/T") ~ 1,
                               c("T/C") ~ 2,
                               c("C/T") ~ 3,
                               c("C/C") ~ 4,
                               .default = NA)
  )



df$PICALM_rs3851179 <- factor(df$PICALM_rs3851179)
df <- dummy_cols(df, select_columns = ("PICALM_rs3851179"), remove_first_dummy = TRUE )





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

#Analisi bivariata

#NEO_NEU e BDI
df %>% ggplot(aes(x=NEO_NEU, y=BDI)) + geom_point() +
      geom_smooth(method = lm, color="red", fill = "#ac31f6", se=TRUE) +
      theme_tufte()

df %>% ggplot(aes(x=NEO_NEU, y=BDI)) + geom_density2d_filled() +
  theme_tufte()

#NEO_NEU e SES
df %>% ggplot(aes(x=NEO_NEU, y=SES)) + geom_point() +
  geom_smooth(method = lm, color="red", fill = "#ac31f6", se=TRUE) +
  theme_tufte()

df %>% ggplot(aes(x=NEO_NEU, y=SES)) + geom_density2d_filled() +
  theme_tufte()

#MINI COPE 7 e 8
df %>% ggplot(aes(x=MINI.COPE_7, y=MINI.COPE_8)) + geom_point() +
  geom_smooth(method = lm, color="red", fill = "#ac31f6", se=TRUE) +
  theme_tufte()

df %>% ggplot(aes(x=MINI.COPE_7, y=MINI.COPE_8)) + geom_density2d_filled() +
  theme_tufte()

df_pca <- df %>%
  select(leukocytes : HSV_r)
df_pca <- na.omit(df_pca)
pca_res <- prcomp(df_pca, scale. = TRUE)

fviz_eig(
  pca_res,
  addlabels = TRUE,
  ncp = 15,
  ylim = c(0, 80),
  main = "Scree plot con percentuali"
)

fviz_contrib(pca_res,
             choice = "var",
             axes = 1:4,         # Somma l'importanza su PC1 & PC2 & PC3 & PC4;
             top = 28)


#Estrazione degli autovalori e delle percentuali
eig_val <- get_eigenvalue(pca_res)

#Creazione di un dataframe leggibile
tabella_varianza <- data.frame(
  Dimensione = 1:nrow(eig_val),
  Varianza_Percentuale = round(eig_val$variance.percent, 2),
  Varianza_Cumulata = round(eig_val$cumulative.variance.percent, 2)
)

# Trova il numero minimo di dimensioni per arrivare al 70%
dim_70 <- which(tabella_varianza$Varianza_Cumulata >= 70)[1]

# Trova il numero minimo di dimensioni per arrivare all'80%
dim_80 <- which(tabella_varianza$Varianza_Cumulata >= 80)[1]

cat("Per spiegare almeno il 70% dell'informazione servono:", dim_70, "dimensioni\n")
cat("Per spiegare almeno il 80% dell'informazione servono:", dim_80, "dimensioni\n")


# Guarda quali variabili pesano di più sulle prime 7 dimensioni
res_var <- get_pca_var(pca_res)
contributi <- res_var$contrib[, 1:7]


sort(contributi[, 1], decreasing = TRUE)[1:5]   # top 5 variabili della PC1:
sort(contributi[, 2], decreasing = TRUE)[1:5]   # top 5 variabili della PC2:
sort(contributi[, 3], decreasing = TRUE)[1:5]   # top 5 variabili della PC3:
sort(contributi[, 4], decreasing = TRUE)[1:5]   # top 5 variabili della PC4:
sort(contributi[, 5], decreasing = TRUE)[1:5]   # top 5 variabili della PC5:
sort(contributi[, 6], decreasing = TRUE)[1:5]   # top 5 variabili della PC6:
sort(contributi[, 7], decreasing = TRUE)[1:5]   # top 5 variabili della PC7:

analysis_name <- df %>% select(leukocytes:HSV_r) %>% colnames()
dataframe_pca <- df %>% select(-all_of(analysis_name))
dataframe_pca <- dataframe_pca[rownames(pca_res$x), ]

# Otteniamo le 7 dimensioni scelte dallo Scree Plot
nuove_dimensioni <- as.data.frame(pca_res$x[, 1:7])
colnames(nuove_dimensioni) <- paste0("Dim_", 1:7)

# Aggiungiamo le nuovi dimensioni ottenute
dataframe_pca <- cbind(dataframe_pca, nuove_dimensioni)

dataframe_pca %>% ggplot( aes(x=Dim_1)) +
  geom_density(fill="#696089", color="#000000", alpha=0.8) +
  ggtitle("Dim 1") + theme_minimal()


ggplot(df, aes(APOE_risk_score, CVLT_13)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.line = element_line(linewidth = 0.75),
    axis.line.x.bottom = element_line(colour = "blue")
  )
