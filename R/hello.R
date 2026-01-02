library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(readr)
library(janitor)
library(corrplot)
library(fastDummies)
library(skimr)
library(psych)
library(factoextra)

dt <- read_tsv(file = file.path(getwd(), "dataset", "participants.tsv"), na = c("", "NA", "n/a", "N/A"))

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

#imputazione istruzione
df <- df %>%
  mutate(education = ifelse(is.na(education), min(education, na.rm = TRUE), education))

#imputazione allergie
df <- df %>%
  mutate(allergies = ifelse(is.na(allergies), min(allergies, na.rm = TRUE), allergies))

#imputazione ibuprofene
df <- df %>%
  mutate(ibuprofen_intake = replace_na(ibuprofen_intake, 0))

df <- df %>%
  mutate(thyroid_diseases = ifelse(is.na(thyroid_diseases), round(median(thyroid_diseases, na.rm = TRUE)),thyroid_diseases))

df <- df %>%
  mutate(smoking_status = ifelse(is.na(smoking_status), round(median(smoking_status, na.rm = TRUE)), smoking_status))

df <- df %>%
  mutate(hypertension = ifelse(is.na(hypertension), round(median(hypertension, na.rm = TRUE)), hypertension))

df <- df %>%
  mutate(other_diseases = ifelse(is.na(other_diseases),round(median(other_diseases, na.rm = TRUE)), other_diseases))


# --- GRAFICO 1: APOE Haplotype ---
ggplot(data = df, aes(x = APOE_haplotype, fill = factor(APOE_haplotype))) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "Genes Presence (APOE Haplotype)",
    x = "Haplotype",
    y = "Count",
    fill = "Genes"
  )

# --- GRAFICO 2: apoe rs429358 ---
ggplot(data = df, aes(x = APOE_rs429358, fill = factor(APOE_rs429358))) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "APOE rs429358 Presence",
    x = "Genotype",
    y = "Count",
    fill = "Genotype"
  )

# --- GRAFICO 3: apoe rs7412 ---
ggplot(data = df, aes(x = APOE_rs7412, fill = factor(APOE_rs7412))) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "APOE rs7412 Presence",
    x = "Genotype",
    y = "Count",
    fill = "Genotype"
  )

# --- GRAFICO 4: Diabete ---
ggplot(data = df, aes(x = diabetes, fill = factor(diabetes))) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "Distribuzione del Diabete",
    x = "Diabetes (0 = No, 1 = Yes)",
    y = "Count",
    fill = "Diabetes"
  )

# --- GRAFICO 5: Ipertensione ---
ggplot(data = df, aes(x = hypertension, fill = factor(hypertension))) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "Distribuzione Hypertension",
    x = "Hypertension (0 = No, 1 = Yes)",
    y = "Count",
    fill = "Hypertension"
  )

# --- GRAFICO 6: BDI Boxplot ---
ggplot(data = df, aes(x = "", y = BDI)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "BDI Boxplot (Beck Depression Inventory)",
    x = "",
    y = "Score"
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

# Creazione del dataset lungo per il grafico
df_long <- pivot_longer(df,
                        cols = c("NEO_NEU", "NEO_EXT", "NEO_OPE", "NEO_AGR", "NEO_CON"),
                        names_to = "Variabile",
                        values_to = "Valore")

# Grafico corretto (senza caratteri nascosti)
ggplot(df_long, aes(x = Valore, fill = Variabile)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(~ Variabile, ncol = 2) +
  theme_minimal() +
  labs(title = "Confronto Densità dei valori NEO", fill = "Variabile")

# --- CORREZIONE BLOCCO RISK SCORE ---
# Ho riscritto questo pezzo rimuovendo gli spazi invisibili che causavano l'errore
df <- df %>%
  mutate(
    APOE_risk_score = case_match(APOE_haplotype,
                                 c("e2/e2", "e3/e2") ~ 1,
                                 "e3/e3"             ~ 2,
                                 c("e3/e4", "e2/e4") ~ 3,
                                 "e4/e4"             ~ 4,
                                 .default = NA)
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





# --- NUOVO CODICE CORRELAZIONE (SENZA ERRORI) ---
library(corrplot)

# 1. Selezioniamo solo le colonne numeriche
df_numeric <- df %>% select(where(is.numeric))

# 2. Rimuoviamo colonne con varianza zero (cioè che hanno sempre lo stesso valore)
#    e gestiamo i valori mancanti mettendo 0
df_clean <- df_numeric[, sapply(df_numeric, var, na.rm = TRUE) > 0]
df_clean[is.na(df_clean)] <- 0

# 3. Selezioniamo solo le prime 20 variabili per rendere il grafico leggibile
df_corr_subset <- df_clean %>% select(1:min(20, ncol(df_clean)))

# 4. Generiamo il grafico stabile
corrplot(cor(df_corr_subset),
         method = "circle",
         type = "upper",
         tl.cex = 0.6,
         title = "Matrice di Correlazione")

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

print(fviz_contrib(pca_res,
                   choice = "var",
                   axes = 1:4,
                   top = 28,
                   title = "Contributo variabili alle prime 4 PC"))

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


# --- VISUALIZZAZIONE CONTRIBUTI VARIABILI (Codice Pulito) ---

cat("\n--- Top 5 variabili per PC1 ---\n")
print(sort(contributi[, 1], decreasing = TRUE)[1:5])

cat("\n--- Top 5 variabili per PC2 ---\n")
print(sort(contributi[, 2], decreasing = TRUE)[1:5])

cat("\n--- Top 5 variabili per PC3 ---\n")
print(sort(contributi[, 3], decreasing = TRUE)[1:5])

cat("\n--- Top 5 variabili per PC4 ---\n")
print(sort(contributi[, 4], decreasing = TRUE)[1:5])

cat("\n--- Top 5 variabili per PC5 ---\n")
print(sort(contributi[, 5], decreasing = TRUE)[1:5])

cat("\n--- Top 5 variabili per PC6 ---\n")
print(sort(contributi[, 6], decreasing = TRUE)[1:5])

cat("\n--- Top 5 variabili per PC7 ---\n")
print(sort(contributi[, 7], decreasing = TRUE)[1:5])

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

print("tutto ok")
