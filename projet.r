library(ade4)


env= read.csv('env.csv', header = TRUE,  sep = ';',  stringsAsFactors = FALSE, row.names = 1)
sp= read.csv('sp.csv', header = TRUE,  sep = ';',  stringsAsFactors = FALSE, row.names = 1)
traits= read.csv('traits.csv', header = TRUE,  sep = ';',  stringsAsFactors = FALSE, row.names = 1)

env 
sp 
traits 

str(env)

    ## Preparation des données

    # Environnement

env=env[,-c(1,2,3,8)]

    # Traits

breaks_troph <- quantile(traits$TROPH, probs = seq(0, 1, 0.25), na.rm = TRUE)

traits$TROPH_classe <- cut(traits$TROPH, 
                      breaks = breaks_troph, 
                      include.lowest = TRUE,
                      labels = c("Q1", "Q2", "Q3", "Q4"))


breaks_TL <- quantile(traits$TL, probs = seq(0, 1, 0.25), na.rm = TRUE)

traits$TL_classe <- cut(traits$TL, 
                      breaks = breaks_TL, 
                      include.lowest = TRUE,
                      labels = c("Q1", "Q2", "Q3", "Q4"))


breaks_TMAX <- quantile(traits$TMAX, probs = seq(0, 1, 0.25), na.rm = TRUE)

traits$TMAX_classe <- cut(traits$TMAX, 
                      breaks = breaks_TMAX, 
                      include.lowest = TRUE,
                      labels = c("Q1", "Q2", "Q3", "Q4"))

traits=traits[,-c(1,2,3,7,11,16)]
head(traits)

traits <- as.data.frame(lapply(traits, as.factor))




    ## Analyse 

# Analyse des Espèces 
afc_esp <- dudi.coa(sp, scannf = FALSE, nf = 2)
scatter(afc_esp, posi = "topright")
title("AFC des Espèces")
inertia.dudi(afc_esp)

# Analyse de l'Environnement 
acp1=dudi.pca(env,center = TRUE, scale = TRUE, scannf = F, nf=2)
barplot(acp1$eig, main="Valeurs propres") 
sum(acp1$eig[1:2]) / sum(acp1$eig) * 100
s.corcircle(acp1$co, clabel = 0.8) 
plot(acp1$li)



# Analyse des traits

acm_traits <- dudi.acm(traits, scannf = FALSE, nf = 2)
s.arrow(acm_traits$co, clabel = 0.7, sub = "Traits fonctionnels (ACM)")





# Analyse conjointe des espèces et de l'environnement

acp_poid=dudi.pca(env,row.w=afc_esp$lw, scannf = FALSE, nf=2)

coiner1=coinertia(acp_poid, afc_esp, scannf = FALSE, nf = 2)
test_coiner=randtest(coiner1, nrepet = 999)
plot(coiner1, main = "Analyse conjointe des espèces et de l'environnement")
print(test_coiner)



# Analyse conjointe des espèces et des traits
