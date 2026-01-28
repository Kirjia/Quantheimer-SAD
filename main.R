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
library(heatmaply)
library(ClusterR)
library(htmlwidgets)

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

#imputazione CVLT_12
df <- df %>%
  mutate(CVLT_12 = replace(CVLT_12, values =median(df$CVLT_12, na.rm = TRUE)))

#imputazione CVLT_13
df <- df %>%
  mutate(CVLT_13 = replace(CVLT_13, values= median(df$CVLT_13, na.rm = TRUE)))

df <- df %>%
  mutate(learning_deficits = replace(learning_deficits, values= round(mean(learning_deficits, na.rm = TRUE))))


df <- df %>%
  mutate(thyroid_diseases = ifelse(is.na(thyroid_diseases), round(median(thyroid_diseases, na.rm = TRUE)),thyroid_diseases))

df <- df %>%
  mutate(smoking_status = ifelse(is.na(smoking_status), round(median(smoking_status, na.rm = TRUE)), smoking_status))

df <- df %>%
  mutate(hypertension = ifelse(is.na(hypertension), round(median(hypertension, na.rm = TRUE)), hypertension))

df <- df %>%
  mutate(other_diseases = ifelse(is.na(other_diseases),round(median(other_diseases, na.rm = TRUE)), other_diseases))

#imputazione MINI.cope
colonne_MINI <- c("MINI.COPE_1", "MINI.COPE_2", "MINI.COPE_3", "MINI.COPE_4", "MINI.COPE_5", "MINI.COPE_6","MINI.COPE_7", "MINI.COPE_8", "MINI.COPE_9","MINI.COPE_10", "MINI.COPE_11", "MINI.COPE_12","MINI.COPE_13", "MINI.COPE_14")

df$participant_id <- as.numeric(gsub("sub-", "", df$participant_id))

df <- df %>%
  mutate(across(all_of(colonne_MINI),
                ~ ifelse(is.na(.), round(median(., na.rm = TRUE)), .)))



df <- df %>%
  mutate(
    # 1. PROBLEM-FOCUSED / ACTIVE (Coping Adattivo)
    # Agire, pianificare e cercare consigli pratici
    Cope_Active = rowSums(select(.,
                                 MINI.COPE_1,  # Active Coping
                                 MINI.COPE_2,  # Planning
                                 MINI.COPE_8   # Instrumental Support
    ), na.rm = TRUE),

    # 2. EMOTION-FOCUSED (Coping Emotivo/Neutro)
    # Gestire le emozioni, cercare conforto, ridere, pregare
    Cope_Emotional = rowSums(select(.,
                                    MINI.COPE_3,  # Positive Reframing
                                    MINI.COPE_4,  # Acceptance
                                    MINI.COPE_5,  # Humor
                                    MINI.COPE_6,  # Religion
                                    MINI.COPE_7,  # Emotional Support
                                    MINI.COPE_11  # Venting
    ), na.rm = TRUE),

    # 3. AVOIDANT / DYSFUNCTIONAL (Coping Disfunzionale)
    # Evitare il problema, bere, negare, incolparsi
    # Questo gruppo correla spesso con declino cognitivo e depressione
    Cope_Avoidant = rowSums(select(.,
                                   MINI.COPE_9,  # Self-Distraction
                                   MINI.COPE_10, # Denial
                                   MINI.COPE_12, # Substance Use
                                   MINI.COPE_13, # Behavioral Disengagement
                                   MINI.COPE_14  # Self-Blame
    ), na.rm = TRUE)
  )

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
df %>% ggplot( aes(x=smoking_status)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione fumatori") + theme_minimal()
df %>% ggplot( aes(x=drugs)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione farmaci assunti") + theme_minimal()
df %>% ggplot( aes(x=ibuprofen_intake)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione ibuprofene") + theme_minimal()
df %>% ggplot( aes(x=coffee_status)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione assunzione caffe") + theme_minimal()
df %>% ggplot( aes(x=thyroid_diseases)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione malattie tiroidee") + theme_minimal()
df %>% ggplot( aes(x=other_diseases)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione altre malattie presenti)") + theme_minimal()
df %>% ggplot( aes(x=Cope_Active)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione code active") + theme_minimal()
df %>% ggplot( aes(x=Cope_Emotional)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione cope emotional") + theme_minimal()
df %>% ggplot( aes(x=Cope_Avoidant)) +
  geom_density(fill="#696089", color="#000000", alpha=0.9) +
  ggtitle("Distribuzione cope avoidant") + theme_minimal()


# Creazione del dataset lungo per il grafico
df_long <- pivot_longer(df,
                        cols = c("NEO_NEU", "NEO_EXT", "NEO_OPE", "NEO_AGR", "NEO_CON"),
                        names_to = "Variabile",
                        values_to = "Valore")


ggplot(df_long, aes(x = Valore, fill = Variabile)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(~ Variabile, ncol = 2) +
  theme_minimal() +
  labs(title = "Confronto Densità dei valori NEO", fill = "Variabile")


select(df, 'CVLT_1':'CVLT_13') %>%
  pivot_longer(cols= 'CVLT_1':'CVLT_13', names_to = "Variabile", values_to = "Valore") %>%
  ggplot(aes(x = Valore, fill = Variabile)) +
  geom_density(alpha = 0.6, color = NA, trim=TRUE) +
  facet_wrap(~ Variabile, ncol = 3, scales = "free") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Confronto Densità CVLT")

select(df, 'CVLT_1':'CVLT_13') %>%
  pivot_longer(cols= 'CVLT_1':'CVLT_13', names_to = "Variabile", values_to = "Valore") %>%
  ggplot(aes(x = Valore, fill = Variabile)) +

  geom_histogram(bins = 15, color = "white", alpha = 0.8) +
  facet_wrap(~ Variabile, ncol = 3, scales = "free") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Distribuzione Punteggi CVLT (Istogramma)")


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



# Verifica rapida
print("Nuove variabili create:")
head(select(df, Cope_Active, Cope_Emotional, Cope_Avoidant))

df <- select(df, -"MINI.COPE_1":-"MINI.COPE_14")


#Grafico della correlazione
df_clean <- select(df, "APOE_rs429358", "APOE_rs7412", "age":"Cope_Avoidant")
df_clean <- df_clean[, apply(df_clean, 2, var, na.rm = TRUE) != 0]
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

heatmap_iniziale <- heatmaply_cor(
  cor(df_clean, use = "pairwise.complete.obs"),
  # Opzionale: Estetica del testo
  node_type = "scatter",
  hclust_method= "ward.D2",
  point_size_mat = -log10(p),
  point_size_name = "-log10(p-value)",
  dendrogram = "none",
  label_names = c("x", "y", "Correlation")
)

print(heatmap_iniziale)

saveWidget(heatmap_iniziale, file="Heatmap iniziale.html", selfcontained = TRUE)

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


ggplot(df, aes(APOE_risk_score, CVLT_9)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.line = element_line(linewidth = 0.75),
    axis.line.x.bottom = element_line(colour = "blue")
  )

ggplot(df, aes(APOE_risk_score, CVLT_7)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.line = element_line(linewidth = 0.75),
    axis.line.x.bottom = element_line(colour = "blue")
  )


# ANALISI QUANTILI#

colonne_numeriche <- df %>% select("age":"APOE_risk_score") %>%
  select(where(is.numeric)) %>%
  names()

func_quantili <- function(df_input, col_name_input_input){


  dati_vector <- df_input[[col_name_input_input]]

  q_vals <- quantile(dati_vector, probs = c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)


  p <- ggplot(df_input, aes(x = .data[[col_name_input_input]])) +
    geom_density(fill = "#69b3a2", alpha = 0.6) +
    geom_rug(alpha = 0.1) +
    geom_vline(xintercept = q_vals[1], color = "red", linetype = "dotted", linewidth=1) +
    geom_vline(xintercept = q_vals[5], color = "red", linetype = "dotted", linewidth=1) +
    geom_vline(xintercept = q_vals[2], color = "blue", linetype = "dashed") +
    geom_vline(xintercept = q_vals[4], color = "blue", linetype = "dashed") +
    geom_vline(xintercept = q_vals[3], color = "black", linewidth = 1.2) +
    labs(
      title = paste("Analisi Quantili:", col_name_input_input),
      subtitle = paste0("Mediana: ", round(q_vals[3], 2)),
      x = col_name_input_input,
      y = "Densità"
    ) +
    theme_minimal()

  print(p)
}

func_boxplot <- function(df_input, col_name_input){

  # A. Calcoliamo gli outlier
  valori <- df_input[[col_name_input]]
  n_outliers <- length(boxplot.stats(valori)$out)

  # B. Creiamo il Boxplot

  p <- ggplot(df_input, aes(x = "", y = .data[[col_name_input]])) +

    # 1. Boxplot
    geom_boxplot(fill = "orange", alpha = 0.6,
                 outlier.colour = "red", outlier.shape = 16, outlier.size = 3) +

    # 2. Jitter
    geom_jitter(width = 0.2, alpha = 0.1, color = "black") +


    theme_minimal() +
    theme(
      axis.title.x = element_blank(), # Toglie il titolo asse X
      axis.text.x = element_blank(),  # Toglie le etichette asse X
      axis.ticks.x = element_blank()
    ) +
    labs(
      title = paste("Boxplot:", col_name_input),
      subtitle = paste("Numero di Outlier rilevati:", n_outliers),
      y = col_name_input
    )

  print(p)
}

report_pdf_safe <- function(colonne_numeriche, nome_file, do_func) {


  if(file.exists(nome_file)) {
    message("Nota: Il file esistente verrà sovrascritto.")
  }

  pdf(nome_file, width = 10, height = 6)

  # ---------------------------------------------------------
  #
  # Questo assicura che il dispositivo grafico venga chiuso
  # appena la funzione finisce, anche se c'è un errore.
  # ---------------------------------------------------------
  on.exit({
    dev.off()
    message("Dispositivo grafico chiuso correttamente. Controllo schermo ripristinato.")
  }, add = TRUE)



  if(length(colonne_numeriche) == 0) {
    stop("Errore: Nessuna colonna numerica trovata nel dataframe.")
  }

  message(paste("Inizio generazione grafici per", length(colonne_numeriche), "colonne..."))


  for (col in colonne_numeriche) {

    try({
      # Controllo rapido se la colonna è vuota prima di chiamare la funzione
      valori <- df[[col]]
      if (all(is.na(valori))) next


      do_func(df, col)

    }, silent = FALSE)
  }

  message(paste("Finito! Il file è stato salvato come:", nome_file))
}
#fine quantili#

report_pdf_safe(colonne_numeriche, "Report Quantili.pdf", func_quantili)
report_pdf_safe(colonne_numeriche, "Report Boxplot.pdf", func_boxplot)

############################################## PCA--------------CVLT----------------PCA #############################################

#Visulizzazione varizanza CVLT tranne 9 per tenere il suo valore separato, 12 e 13 in pratica costanti quindi verranno rimosse
CVLT_pca <- select(df, 'CVLT_1':'CVLT_11') %>% select(-'CVLT_9') %>%
  prcomp(scale. = TRUE)
fviz_eig(CVLT_pca, addlabels = TRUE, ncp = 13, ylim=c(0, 80), main="CVLT Var distr.")
print(fviz_contrib(CVLT_pca,
                   choice = "var",
                   axes = 1:4,
                   top = 28,
                   title = "Contributo variabili alle prime 4 PC"))

#Estrazione degli autovalori e delle percentuali
CVLT_eig_val <- get_eigenvalue(CVLT_pca)

#Creazione di un dataframe leggibile
CVLT_tabella_varianza <- data.frame(
  Dimensione = 1:nrow(CVLT_eig_val),
  Varianza_Percentuale = round(CVLT_eig_val$variance.percent, 2),
  Varianza_Cumulata = round(CVLT_eig_val$cumulative.variance.percent, 2)
)
# Otteniamo le dimensioni scelte dallo Scree Plot
nuove_dimensioni <- as.data.frame(CVLT_pca$x[, 1:2])
colnames(nuove_dimensioni) <- paste0("N_CVLT_", 1:2)

# Aggiungiamo le nuovi dimensioni ottenute
df <- cbind(df, nuove_dimensioni)

df <- select(df, -c("CVLT_10", "CVLT_11", "CVLT_12", "CVLT_13", "learning_deficits", "diabetes", "APOE_rs429358", "APOE_rs7412"))
df <- select(df, -"CVLT_1":-"CVLT_8")



############################################## PCA--------------END CVLT----------------PCA #############################################

############################################## PCA--------------Blood Test----------PCA #############################################
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


# --- VISUALIZZAZIONE CONTRIBUTI VARIABILI ---

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

fviz_contrib(pca_res, choice = "var", axes = 1, top = 10)

# Otteniamo le 5 dimensioni scelte dallo Scree Plot
nuove_dimensioni <- as.data.frame(pca_res$x[, 1:5])
colnames(nuove_dimensioni) <- paste0("Dim_", 1:5)

# Aggiungiamo le nuovi dimensioni ottenute
dataframe_pca <- cbind(dataframe_pca, nuove_dimensioni)

dataframe_pca %>% ggplot( aes(x=Dim_1)) +
  geom_density(fill="#696089", color="#000000", alpha=0.8) +
  ggtitle("Dim 1") + theme_minimal()






############################################## PCA--------------Blood Test----------PCA #############################################
df <- select(df, -c('APOE_haplotype', 'PICALM_rs3851179', "second_phase", "session_order", "allergies", "hypertension", "education", "other_diseases", "thyroid_diseases", "smoking_status", "drugs","ibuprofen_intake"))
df<- select(df, -'leukocytes':-'HSV_r')
dataframe_pca <- select(dataframe_pca, -c('APOE_haplotype', 'PICALM_rs3851179', "second_phase", "session_order", "allergies", "hypertension", "education"))

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = BMI, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = BMI, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +

  theme_minimal() +
  labs(title = "Confronto BMI", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = APOE_risk_score, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = APOE_risk_score, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +

  theme_minimal() +
  labs(title = "Confronto AOPE risk score", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = BDI, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = BDI, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto BDI", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = SES, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = SES, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto SES", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = RPM, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = RPM, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto RPM", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = EHI, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = EHI, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto EHI", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = NEO_NEU, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = NEO_NEU, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto NEO-NEU", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = NEO_EXT, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = NEO_EXT, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto NEO-EXT", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = NEO_OPE, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = NEO_OPE, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto NEO-OPE", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = NEO_AGR, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = NEO_AGR, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto NEO-AGR", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = NEO_CON, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = NEO_CON, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto NEO-CON", x = "Valore", y = "Densità")


ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = AUDIT, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = AUDIT, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto AUDIT", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = CVLT_9, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = CVLT_9, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto CVLT 9(long free recall delay)", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = sex, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = sex, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto sex", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = age, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = age, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto age", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = dementia_history_parents, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = dementia_history_parents, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto storia demenza familiare", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = `PICALM_rs3851179_G/A`, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = `PICALM_rs3851179_G/A`, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto gene PICALM", x = "Valore", y = "Densità")

ggplot() +
  # Primo dataset
  geom_density(data = df, aes(x = coffee_status, fill = "full samples"), alpha = 0.4) +
  # Secondo dataset
  geom_density(data = dataframe_pca, aes(x = coffee_status, fill = "blood samples"), alpha = 0.4) +
  scale_fill_manual(name = "Legenda",
                    values = c("full samples" = "blue", "blood samples" = "red")) +
  theme_minimal() +
  labs(title = "Confronto coffee status", x = "Valore", y = "Densità")



df <- select(df, -"NEO_NEU", -"NEO_EXT", -"coffee_status", -"PICALM_rs3851179_G/G")


id_col <- "participant_id"

colonne_da_aggiungere <- select(dataframe_pca, "Dim_1":"Dim_5") %>% names() %>% setdiff(names(df))
df_sangue_clean <- dataframe_pca %>%
  select(all_of(id_col), all_of(colonne_da_aggiungere))

dataframe_alz <- inner_join(df, df_sangue_clean, by=id_col)


df_clust <- select(dataframe_alz, -c("participant_id", "sex", "APOE_risk_score", "dementia_history_parents", "PICALM_rs3851179_G/A"))
df_scaled <- scale(df_clust)

p <- cor.test.p(df_scaled)

heatmap_finale <- heatmaply_cor(
  cor(df_scaled, use = "pairwise.complete.obs"),

  node_type = "scatter",
  hclust_method= "ward.D2",
  point_size_mat = -log10(p),
  point_size_name = "-log10(p-value)",
  dendrogram = "none",
  label_names = c("x", "y", "Correlation")
)

print(heatmap_finale)

saveWidget(heatmap_finale, file="Heatmap finale.html", selfcontained = TRUE)

remove(df_long)
remove(df_sangue_clean)
remove(df_clean)
remove(nuove_dimensioni)
remove(dataframe_pca)

###########BEGIN#####CLUSTER#################BEGIN#####CLUSTER###################BEGIN#####CLUSTER#########################
set.seed(123)
km_pp_func <- function(x, k) {
  km <- KMeans_rcpp(x, clusters = k, initializer = 'kmeans++', num_init = 25)
  return(list(cluster = km$clusters, tot.withinss = sum(km$WCSS_sum)))
}


fviz_nbclust(df_scaled,
             FUNcluster = km_pp_func,
             method="wss") +
  geom_vline(xintercept = 3, linetype = 2, color = "red") +
  labs(title = "Metodo Elbow con k-means++",
       subtitle = "Somma dei quadrati entro i cluster per ogni k")


#usando K= 3

km_rcpp_2 <- KMeans_rcpp(df_scaled, clusters = 2, initializer = 'kmeans++', num_init = 25)

fviz_nbclust(df_scaled, FUNcluster = km_pp_func, method = "silhouette")

# 2. Visualizza con fviz_cluster
# Dobbiamo passare i dati e le assegnazioni dei cluster
fviz_cluster(list(data = df_scaled, cluster = km_rcpp_2$clusters),
             ellipse.type = "convex",
             palette = "jco",
             geom = c("point"),
             ggtheme = theme_minimal())

km_rcpp_3 <- KMeans_rcpp(df_scaled, clusters = 3, initializer = 'kmeans++')

fviz_nbclust(df_scaled, FUNcluster = km_pp_func, method = "silhouette")

# 2. Visualizza con fviz_cluster
# Dobbiamo passare i dati e le assegnazioni dei cluster
fviz_cluster(list(data = df_scaled, cluster = km_rcpp_3$clusters),
                         ellipse.type = "convex",
                         palette = "jco",
                         geom = c("point"),
                         ggtheme = theme_minimal())

df_clust <- select(dataframe_alz, c("participant_id", "sex", "APOE_risk_score", "dementia_history_parents", "PICALM_rs3851179_G/A")) %>% cbind(df_clust)
df_clust$Cluster <- as.factor(km_rcpp_3$clusters)

wcss_2 <- km_rcpp_2$WCSS_per_cluster
global_mean <- colMeans(df_scaled)
tss <- sum(apply(df_scaled, 1, function(x) sum((x - global_mean)^2)))
bcss_2 <- tss - wcss_2

n <- nrow(df_scaled)
k <- 2

ch_index_2 <- (bcss_2 / (k - 1)) / (wcss_2 / (n - k))


wcss_3 <- km_rcpp_3$WCSS_per_cluster
bcss_3 <- tss - wcss_3
k <- 3

ch_index_3 <- (bcss_3 / (k - 1)) / (wcss_3 / (n - k))

cat("--- Risultati Clustering ---\n")
cat("WCSS (Within-Cluster SS):  2 Cluster ", wcss_2, " 3 Cluster", wcss_3, "\n")
cat("BCSS (Between-Cluster SS): 2 Cluster ", bcss_2, " 3 cluster ", bcss_3, "\n")
cat("TSS (Total SS):          ", tss, "\n")
cat("Indice Calinski-Harabasz: 2 Cluster ", ch_index_2, " 3 CLuster ", ch_index_3, "\n")

df_summary <- df_clust %>%
  group_by(Cluster) %>%
  summarise(
    N_Pazienti = n(),
    CVLT_1 = mean(N_CVLT_1),
    Hematic_Constitution = mean(Dim_1),
    Depressione_Media = mean(BDI),
    Coping_Evitante = mean(Cope_Avoidant),
    Raven_RPM = mean(RPM),
    Eta = mean(age),

    apoe_risk = mean(APOE_risk_score == 4) * 100
  )

print(df_summary)

df_plot_long <- df_clust %>%

  select(
    Cluster,
    "Memoria (CVLT Dim1)" = N_CVLT_1,
    "Ossigenazione & Immunità (Blood Dim1)" = Dim_1,
    "Metabolismo (Blood Dim2)" = Dim_2,
    "Depressione (BDI)" = BDI,
    "Riserva Cognitiva (RPM)" = RPM,
    "Coping Evitante" = Cope_Avoidant,

  ) %>%

  pivot_longer(
    cols = -Cluster,
    names_to = "Variabile",
    values_to = "Punteggio"
  )


ggplot(df_plot_long, aes(x = Cluster, y = Punteggio, fill = Cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +


  facet_wrap(~ Variabile, scales = "free_y", ncol = 2) +


  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  ) +
  labs(
    title = "Confronto dei Profili Clinici tra i Cluster",
    subtitle = "Distribuzione delle variabili chiave per ogni gruppo identificato",
    y = "Valore Standardizzato (Z-Score)",
    x = "Gruppo (Cluster)"
  ) +
  scale_fill_brewer(palette = "Set2")


df_plot_long_2 <-df_clust %>% select(
  Cluster,
  "APOE risk" = APOE_risk_score,
  "cvlt 9" = CVLT_9,
  "sesso" = sex,
  "storia demenza familiare" = dementia_history_parents
) %>%

  pivot_longer(
    cols = -Cluster,
    names_to = "Variabile",
    values_to = "Punteggio"
  )


ggplot(df_plot_long_2, aes(x = Cluster, y = Punteggio, fill = Cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +


  facet_wrap(~ Variabile, scales = "free_y", ncol = 2) +


  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 11), # Titoli dei riquadri in grassetto
    legend.position = "bottom"
  ) +
  labs(
    title = "Confronto dei Profili Clinici tra i Cluster",
    subtitle = "Distribuzione delle variabili chiave per ogni gruppo identificato",
    y = "Valore Standardizzato (Z-Score)", # Se hai usato dati scalati
    x = "Gruppo (Cluster)"
  ) +
  scale_fill_brewer(palette = "Set2")

table(df_clust$Cluster, df_clust$sex)
table(df_clust$Cluster, df_clust$APOE_risk_score)




centroids <- km_rcpp_3$centroids


colnames(centroids) <- colnames(df_scaled)


centroids_long <- as.data.frame(centroids) %>%
  mutate(Cluster = paste("Cluster", 1:nrow(centroids))) %>%
  pivot_longer(-Cluster, names_to = "Variabile", values_to = "Valore_Z")

# 4. Visualizza: Le barre più alte
ggplot(centroids_long, aes(x = Variabile, y = Valore_Z, fill = Cluster)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Impatto delle variabili per Cluster (Z-Score)",
       y = "Distanza dalla media globale (0)",
       x = "Variabile") +
  theme_minimal()


write.csv2(df_clust, "Quantum/quantum-alz.csv")
###########END#####CLUSTER#################END#####CLUSTER###################END#####CLUSTER#########################





