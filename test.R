# ==============================================================================
# PROJET : Syndrome des cours d'eau & Homogénéisation Biotique (Gagnants/Perdants)
# ==============================================================================

library(ade4)

# ==============================================================================
# ÉTAPE 0 : PRÉPARATION DES DONNÉES
# ==============================================================================
env <- read.csv("env.csv", sep=";", row.names=1)
sp <- read.csv("sp.csv", sep=";", row.names=1)
traits <- read.csv("traits.csv", sep=";", row.names=1)

# Nettoyage des zéros absolus
sp <- sp[, colSums(sp) > 0]
sp <- sp[rowSums(sp) > 0, ]
env <- env[rownames(sp), ]

# Variables environnementales (transformation log sur la morphologie)
env_quant <- env[, c("TEMP", "MTC", "MTW", "PREC", "BV", "DIST", "ALT", "LARG", "PROF")]
env_quant$BV   <- log(env_quant$BV + 1)
env_quant$DIST <- log(env_quant$DIST + 1)
env_quant$ALT  <- log(env_quant$ALT + 1)
env_quant$LARG <- log(env_quant$LARG + 1)
env_quant$PROF <- log(env_quant$PROF + 1)

# Traits ciblés sur la VULNÉRABILITÉ
traits_vuln <- traits[, c("QUAL", "HAB", "OXY", "TEMP", "STATUT")]
traits_vuln <- as.data.frame(lapply(traits_vuln, as.factor))


# ==============================================================================
# ÉTAPE 1 : LE PORTRAIT DES GAGNANTS ET PERDANTS (ACM + CAH)
# ==============================================================================
acm_vuln <- dudi.acm(traits_vuln, scannf = FALSE, nf = 2)

dist_vuln <- dist(acm_vuln$li)
cah_vuln <- hclust(dist_vuln, method = "ward.D2")
traits$Profil <- as.factor(cutree(cah_vuln, k = 3))

par(mfrow = c(1, 2))
s.label(acm_vuln$co, clabel = 0.8)
plot(cah_vuln, labels = rownames(traits), cex = 0.6, xlab = "", sub = "")
rect.hclust(cah_vuln, k = 3, border = c("red", "orange", "green3"))
par(mfrow = c(1, 1))


# ==============================================================================
# ÉTAPE 2 : LE SYNDROME ANTHROPIQUE DU MILIEU (ACP)
# ==============================================================================
pca_env <- dudi.pca(env_quant, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# Projection de l'Occupation du Sol sur l'ACP environnementale
s.class(pca_env$li, fac = as.factor(env$OS), 
        col = c("gold", "red", "forestgreen", "blue"),
        cstar = 1, cellipse = 1.5, axesell = FALSE, clabel = 0.8)


# ==============================================================================
# ÉTAPE 3 : LE TURNOVER SPÉCIFIQUE (AFC)
# ==============================================================================
afc_sp <- dudi.coa(sp, scannf = FALSE, nf = 2)


# ==============================================================================
# ÉTAPE 4 : LA PREUVE DU FILTRE ANTHROPIQUE (CO-INERTIE)
# ==============================================================================
pca_env_pondd <- dudi.pca(env_quant, row.w = afc_sp$lw, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
coia_syndrome <- coinertia(pca_env_pondd, afc_sp, scannf = FALSE, nf = 2)

test_coia <- randtest(coia_syndrome, nrepet = 999)
print("=== TEST DE CO-INERTIE ===")
print(test_coia)

# Visualisation globale des flèches de co-inertie
plot(coia_syndrome)


# ==============================================================================
# ÉTAPE 5 : LES GRAPHIQUES "WAOUH" (Preuve par les Proportions)
# ==============================================================================

# 1. Nettoyage de la modalité "WATB" (Surfaces en eau, trop rares)
filtre_OS <- env$OS != "WATB"
OS_propre <- droplevels(as.factor(env$OS[filtre_OS]))

# 2. Calcul des PROPORTIONS de chaque profil par station
matrice_appartenance <- model.matrix(~ 0 + Profil, data = traits)
colnames(matrice_appartenance) <- paste0("Profil_", levels(traits$Profil))

abondance_profils <- as.data.frame(as.matrix(sp) %*% matrice_appartenance)
proportions_profils <- (abondance_profils / rowSums(abondance_profils)) * 100 

data_visu <- data.frame(OS = OS_propre, proportions_profils[filtre_OS, ])

# --- VISUALISATION 1 : La Co-inertie selon l'Occupation du Sol ---
# L'axe de la co-inertie sépare-t-il les milieux forestiers des milieux agricoles ?
coord_coia <- coia_syndrome$lX[filtre_OS, ]
s.class(coord_coia, fac = OS_propre, 
        col = c("gold", "red", "forestgreen"), 
        cstar = 1, cellipse = 1.5, axesell = FALSE, clabel = 1)


# --- VISUALISATION 2 : Les Boxplots de substitution (Gagnants vs Perdants) ---
# ATTENTION : Remplace "Profil_1" et "Profil_3" par les vrais numéros de tes profils 
# identifiés grâce à la CAH (ex: Profil_1 = Perdants, Profil_3 = Gagnants).

par(mfrow = c(1, 2))

# Boxplot pour les "Perdants" (Espèces sensibles)
boxplot(Profil_1 ~ OS, data = data_visu, 
        col = c("gold", "red", "forestgreen"),
        ylab = "Proportion dans la communauté (%)", 
        xlab = "Occupation du sol",
        las = 1)

# Boxplot pour les "Gagnants" (Espèces opportunistes/exotiques)
boxplot(Profil_3 ~ OS, data = data_visu, 
        col = c("gold", "red", "forestgreen"),
        ylab = "Proportion dans la communauté (%)", 
        xlab = "Occupation du sol",
        las = 1)

par(mfrow = c(1, 1))
# ==============================================================================
# ÉTAPE 6 (BONUS) : LA CARTE GÉOGRAPHIQUE DU SYNDROME
# ==============================================================================

# 1. Extraction des coordonnées spatiales (en gardant le filtre pour enlever WATB)
xy <- env[filtre_OS, c("X", "Y")]

par(mfrow = c(1, 2))

# --- CARTE 1 : L'occupation du sol dans l'espace ---
# On projette les X et Y. cstar=0 et cellipse=0 permettent de n'afficher que les points purs.
par(mfrow=c(1,2))
s.class(xy, fac = OS_propre, 
        col = c("gold", "red", "forestgreen"), 
        cstar = 0, cellipse = 0, axesell = FALSE, grid = FALSE)

# --- CARTE 2 : La répartition géographique des poissons "Gagnants" ---
# La fonction s.value dessine des carrés (positifs) ou cercles (négatifs) 
# dont la taille est proportionnelle au pourcentage de la guilde.
# ATTENTION : Remplace "Profil_3" par le bon numéro de tes Gagnants !
s.value(xy, data_visu$Profil_3, 
        grid = FALSE,
        csize = 0.5) # csize règle la taille des carrés

par(mfrow = c(1, 1))
q


