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
classif_indic <- hclust(dist_indic, method = "ward.D2") # à verifier avec d'autres méthodes de classification hiérarchique
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
