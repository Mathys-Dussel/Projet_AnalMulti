# Comment l'anthropisation des cours d'eau influence-t-elle les traits et la composition des poissons ?

# H1: Les espèces généralistes ont plus de succés que les spécialistes
# H2: Les espèces gagnantes ont des strategies r.
# H3: 

library(ade4)


##################### Préparation des données ##################
env= read.csv('env.csv', header = TRUE,  sep = ';',  stringsAsFactors = FALSE, row.names = 1)
sp= read.csv('sp.csv', header = TRUE,  sep = ';',  stringsAsFactors = FALSE, row.names = 1)
traits= read.csv('traits.csv', header = TRUE,  sep = ';',  stringsAsFactors = FALSE, row.names = 1)

env 
sp 
traits 

env <- env[rownames(sp), ]

str(traits)
traits_indic=  traits[, c("FE","QUAL", "HAB", "OXY", "TEMP", "STATUT")]
traits_indic <- as.data.frame(lapply(traits_indic, as.factor))


env_quant <- env[, c("TEMP", "MTC", "MTW", "PREC", "BV", "DIST", "ALT", "LARG", "PROF")]
env_quant$BV   <- log(env_quant$BV + 1)
env_quant$DIST <- log(env_quant$DIST + 1)
env_quant$ALT  <- log(env_quant$ALT + 1)
env_quant$LARG <- log(env_quant$LARG + 1)
env_quant$PROF <- log(env_quant$PROF + 1)


sp = log(sp + 1)

##################### Ségrégration des espèces  ##################


acm_indic <- dudi.acm(traits_indic, scannf = FALSE, nf = 2)
s.arrow(acm_indic$co, clabel = 0.8)  

dist_indic <- dist(acm_indic$li)
classif_indic <- hclust(dist_indic, method = "average") # à verifier avec d'autres méthodes de classification hiérarchique
traits$Profil <- as.factor(cutree(classif_indic, k = 3)) # 3 pour gagnants / perdants / intermédiaires

plot(classif_indic, labels = rownames(traits), cex = 0.6)
rect.hclust(classif_indic, k = 3, border = c("red", "orange", "forestgreen"))
       

######################### ACP environnementale  ##################

acp_env = dudi.pca(env_quant, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
sum(acp_env$eig[1:2]/sum(acp_env$eig)*100) 

par(mfrow = c(1, 2))

s.corcircle(acp_env$co, xax = 1, yax = 2, clabel = 0.8)

s.class(acp_env$li, fac = as.factor(env$OS), 
        col = c("gold", "red", "forestgreen", "blue"),cstar = 1,cellipse = 1.5, axesell = FALSE, clabel = 1)

######################## AFC avec les espèces  ##################

afc_sp = dudi.coa(sp, scannf = FALSE, nf = 2)

scatter(afc_sp)
inertia.dudi(afc_sp)

########################## Coinertie entre les espèces et les variables environnementales  ##################

acp_pond = dudi.pca(env_quant, row.w = afc_sp$lw, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
coiner_syndorme= coinertia(acp_pond, afc_sp, scannf = FALSE, nf = 2)

randtest(coiner_syndorme, nrepet = 999)

plot(coiner_syndorme)

##################### Réponse à H1 : Les espèces généralistes ont plus de succès que les spécialistes  ##################

# on dégage WATB et tout en facteur
filtre_OS <- env$OS != "WATB"
OS_propre <- droplevels(as.factor(env$OS[filtre_OS]))

# On fait une matrice profil en fct d'esp
matrice_sp_os <- model.matrix(~ 0 + Profil, data = traits)

#matric entre esp et os
profils_os <- as.data.frame(as.matrix(sp) %*% matrice_sp_os)

#  tout en porcentage
prop_profils <- (profils_os / rowSums(profils_os)) * 100 

# le df utiles pour les plot
data_boxplot <- data.frame(OS = OS_propre, prop_profils[filtre_OS, ])



# on refait un plot de coiner sans WATB
coord_coia <- coiner_syndorme$lX[filtre_OS, ]
s.class(coord_coia, fac = OS_propre, 
        col = c("gold", "red", "forestgreen"), 
        cstar = 1, cellipse = 1.5, axesell = FALSE, clabel = 1)



# les plots de boxplot pour les profils

par(mfrow = c(1,3))
boxplot(Profil1 ~ OS, data = data_boxplot, 
        col = c("gold", "red", "forestgreen"),
        ylab = "Proportion dans la communauté (%)", 
        xlab = "Occupation du sol profil 1",
        las = 1)

boxplot(Profil2 ~ OS, data = data_boxplot, 
        col = c("gold", "red", "forestgreen"),
        ylab = "Proportion dans la communauté (%)", 
        xlab = "Occupation du sol profil 2",
        las = 1)

boxplot(Profil3 ~ OS, data = data_boxplot, 
        col = c("gold", "red", "forestgreen"),
        ylab = "Proportion dans la communauté (%)", 
        xlab = "Occupation du sol profil 3",
        las = 1)
par(mfrow = c(1,1))


############### la carto pcq c'est joliiiiii ##########

xy <- env[filtre_OS, c("X", "Y")]

s.class(xy, fac = OS_propre, 
        col = c("gold", "red", "forestgreen"),
        cstar = 0, cellipse = 0, axesell = FALSE, grid = FALSE)


par(mfrow = c(1, 3))
s.value(xy, data_boxplot$Profil1, 
        grid = FALSE,
        csize = 0.5) 

s.value(xy, data_boxplot$Profil2, 
        grid = FALSE,
        csize = 0.5) 
    
s.value(xy, data_boxplot$Profil3, 
        grid = FALSE,
        csize = 0.5) 

par(mfrow = c(1, 1))


library(ggplot2)

# Création de la carte épurée
ggplot(data.frame(x = xy$X, y = xy$Y, OS = OS_propre), aes(x = x, y = y, color = OS)) +
  geom_point(size = 3, alpha = 0.8) + # Points légèrement plus gros pour compenser l'absence d'axes
  scale_color_manual(values = c("gold", "red", "forestgreen")) + 
  coord_fixed() + # Très important pour ne pas déformer la géographie
  theme_void() +  # Supprime TOUT (axes, fond, grille, texte)
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 10))
  ) +
  labs(title = "Structure spatiale de l'occupation du sol",
       color = "Type d'OS")

