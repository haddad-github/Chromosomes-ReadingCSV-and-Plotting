import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_data(file_name, chromosome):
    
    data = pd.read_csv(file_name, sep="\t", usecols=["chromosome", "position", "low_confidence_variant", "pval"]) #lire fichier avec colonnes specifiques
    
    data["pval"] = data["pval"].replace(0, 1e-274) #remplace les 0 dans pval par 1e-274 pour eviter division par 0
    
    data["log_10"] = -np.log10(data["pval"]) #creer nouvelle colone basee sur un calcul d'une autre rangee
    
    data = data.dropna() #enlever toutes les colones vides
    
    data = data.drop(data.index[data["low_confidence_variant"] == 'false']) #BONUS
    
    return data.loc[data["chromosome"] == chromosome] #definir chromosome de l'appel de la fonction

femelle = read_data("vitamin_d.females.tsv.gz", 1) #lire femelle
male = read_data("vitamin_d.males.tsv.gz", 1) #lire male

plt.figure(figsize=(12,4)) #definir dimensions de la figure
plt.subplots_adjust(hspace = 0) #eliminer l'espace vertical

plt.subplot(211) #pour separer les tableaux en 2
plt.scatter(male["position"], male["log_10"], s=1, c="blue", label = "Males") #x, y, taille points, couleur points, etiquette
plt.title("Vitamin D (nmol/L)", weight = "bold") #titre
plt.xticks([]) #enlever chiffres de l'axe des x
plt.ylabel(r"Males $-log_{10}(pval)$", loc = "top") #mettre haut pour faire place a femelle aussi 
plt.axhline(-np.log10(5e-8), color="darkgray", linestyle="--") #ligne de style -- a position indiquee

log10negatif = -femelle["log_10"] #valeurs negative de la colonne log_10

plt.subplot(212)
plt.scatter(femelle["position"], log10negatif, s=1, c="blue", label = "Females")
plt.xlabel("Chromosome 1 positions") #etiquette axe x
plt.ylabel(r"Females $log_{10}(pval)$", loc = "bottom") #etiquette axe y
plt.axhline(y=np.log10(5e-8), color="darkgray", linestyle="--")

plt.savefig("miami.png", dpi=300, bbox_inches="tight")

#plt.show() #pour tester
