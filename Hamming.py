#-*- encoding:utf-8 -*-
#=========================================================================================
# Importation des bibliothèque "nupmy" pour les calculs scientifique et "random" pour de choix
# aléatoire
#=========================================================================================
import numpy as np
import random as ra
#=========================================================================================
#
#
#
#                               PARTIE I: CHIFFREMENT
#
#
#
#=========================================================================================
#       On définit une fonction "conversion binaire" pour convertir en binaire
# elle prend en entrée la chaine de caractère à convertir
#=========================================================================================
def conversion_binaire(x):
#  On crée une chaine vide: "binaire" 
     binaire=''
     for i in x:
          # la fonction "np.binary_repr" de numpy permet la conversion en binaire
          # la fonction "ord" permet la recuperation du code ascii d'une lettre
          # "binaire+=" permet de remplir la variable binaire 
          binaire+=np.binary_repr(ord(i), 8)
     return binaire 
#=========================================================================================
#    On définit une fonction "découpage" pour découper ce qu'on a converti en binaire
# elle prend en entrée la chaine de caractere à découper "x" et la longueur de decoupage qui est "taille"
# "x" doit être une chaine de bits qu'on vient d'obtenir à partir de la fonction "conversion_binaire"
#=========================================================================================
def decoupage(x):
     # on récupère la taille de la variable à découper
     x1=len(x)
     # On crée une liste vide "v=[]"
     v=[]
     # On crée une variable "x2" et on l'initialise à 0
     x2=0
     while x2 < x1:
         # cette ligne nous permet le decoupage de "x" en bloc de longueur égale à "taille"
         x3=x[x2 : x2+4]
         # création d'une liste de string constituée des éléments de x3
         x4=list(x3)
         # "map(int, x4)" permet de transformer cette de string en une liste des int
         x5=map(int, x4)
         # "v.append" permet de remplir la liste "v" qu'on a crée avant la boucle while
         v.append(x5)
         x2+=4
     return v
#=========================================================================================
# On définit une fonction "systematisation" permettant de systématiser la matrice de contrôle 
# NB: cette fonction retourne juste la partie de la matrice différente de la matrice identité
#=========================================================================================
def systematisation(A):
     # "np.shape" permet de récupérer la taille de la matrice
     (t1, t2)=np.shape(A)
     # "t3" n'est rien d'autre que la dimension de notre matrice generatrice obtenue en calculant
     # en réalité: 2**r - r - 1
     t3=t2 - t1
     # On met la matrice A dans une variable pour le manipuler
     mat=A
     # On initialise une matrice rempli de zéro de taille t1xt3
     # qu'on y mettera la partie de la matrice de systématisée différente de la matrice identité
     mat1=np.zeros([t1, t3], dtype=int)
     # On crée une variable temporaire "tmp1" qui permet la permutation de ligne de notre 
     # matrice pendant la systematisation si le besoin se presente
     tmp1=np.ones([1, t2], dtype=int)
     # Cette boucle nous permet de parcourir les lignes de notre matrice à systématiser  
     for i in xrange(t1):
         # On teste si l'élément de la position (i, t3+i) est nul
         if mat[i,t3 + i]==0:
              # Si oui, on crée une variable test qu'on initialise à zéro
              test=0
              # Cette boucle nous permet de parcourir les éléments de la colonne t3 + i
              # pour chercher une position dont l'élément est non nul
              for l in xrange(t1):
                   # On teste si l'élément de la position (l, t3+i) est non nul 
                   if h_1[l,t3 + i]<>0:
                        # Si oui on change la valeur de test et on arrete la boucle
                        test=1
                        break
              # Si test=1 on permute les ligne 
              if test==1:
                   # Ici c'est la procédure qui nous permet la permutation de lignes
                   # On met l'une de ligne i dans la variable temporaire qu'on a créé 
                   tmp1[0,:]=mat[i,:]
                   # On change le contenu de la ligne i avec celui de la ligne l
                   mat[i,]=mat[l,]
                   # On reprend dans la ligne l le contenu de la ligne i qu'on a gardé dans tmp1
                   mat[l,]=tmp1
              else:
                   # si non on dispose d'une colonne nulle donc il y'a et on affiche que 
                   # la matrice ne peut pas être systématisée en utilisant la fonction 
                   # "raise TypeError"
                   raise TypeError('Votre matrice ne peut pas être systématisée')
         # cette boucle nous permet de mettre à zéro tout les autre élément de la colonne t3+i
         # On parcourt la colonne
         for k in xrange(t1):
              # On teste si l'élément de la position (k, t3+i) est non nul
              if (mat[k,t3 + i]<>0)and(k!=i):
                   # Si oui, on fait une somme direct avec ligne i pour mettre l'élément 
                   # de la position (k, t3+i) à zéro
                   mat[k,]^=mat[i,]
     # On récupère la partie la partie de la matrice différente de la matrice identité dans "mat1"
     # en parcourant les colonnes
     for i in xrange(t3):
          mat1[:, i]=mat[:,i]
     return mat1
#=========================================================================================
#  On définit une fonction "matrice_de_control_de_hamm" pour créer la matrice de controle
# elle prend en entrée r qui est le paramètre de hamming
#=========================================================================================
def matrice_de_control_de_hamm(r):
     # On calcul la longueur de notre code
     n=2**r-1
     #construction de la matrice de parité initialisée à zéro
     H=np.zeros([r, n], dtype=int)
     # Cette boucle nous permet de remplir notre matrice notre matrice de parité qu'on à 
     # initialisé avec un contenu rempli de zero
     for i in xrange(1, n+1):
          # Les colonnes de la matrice de controle de hamming sont constituées
          # des nombres allant de 1 à n=2**r - 1 qu'on convertit en binaire dont on dispose
          # de manière décroissante
          # Ici, on convertit en binaire i sur r bit et on met dans la variable x1
          x1=np.binary_repr(i, r)
          # après conversion, crée une liste de string constituée des éléments de x1  
          x2=list(x1)
          # On convertit la liste x2 de string en liste de int avec la commande "map(int, x2)"
          x3=map(int, x2)
          # On remplace maintenant les zeros de l colonne i-1 par x3
          H[:, i-1]=x3
     return H
#=========================================================================================
#  On définit une fonction "matrice_de_Gene_de_hamm" pour créer la matrice generatrice 
#=========================================================================================
def matrice_de_Gene_de_hamm(r):
     # On calcule la longueur
     n=2**r-1
     # On calcule la dimension
     t3=n-r
     # On calcule la matrice de controle en utilisant la fonction "matrice_de_control_de_hamm"
     H=matrice_de_control_de_hamm(r)
     # On sytématise la matrice de controle en utilisant la fonction "systematisation"
     G1=systematisation(H)
     # Ces deux dernière permettent de concatener la matrice identité d'ordre (k,k) et la 
     # transposée de la partie de la matrice de controle qui est différente de l'identité     
     G=np.concatenate((np.eye(t3, dtype=int), G1))
     G=np.transpose(G)
     return G
#=========================================================================================
#         On définit une fonction codage permettant de coder sans ajouter l'erreur 
# Qui prend en entrée une trois variable "mot" à coder qui doit être une liste et G 
# la matrice generatrice et n le nombre des colonnes
#=========================================================================================
def codage(mot, G):
     # On fait produit mot et G grace à la fonction "np.dot" 
     x2=np.dot(mot, G)
     # après le produit il se peut qu'il ait des élément plus grand que 1 donc il faut 
     # On transforme la matrice obtenue après le produit en une liste 
     x3=[x2[i]%2 for i in xrange(7)]
     return x3
#=========================================================================================
#        On définit une fonction "codage2" permettant de coder et ajouter en même temps
# elle prend en entrée 'x' le caractère à coder, A la matrice generatrice, k la dimension 
# du code n la longueur du code
#=========================================================================================
def codage2(x, A):
    # On met 'x' dans une variable qu'on appelle mot 
    mot=x
    # On convertit mot en binaire en utilisant "conversion_binaire"
    x1=conversion_binaire(mot)
    # On découpe le mot converti en binaire en utilisant "decoupage" qui nous donne une liste de listes
    x2=decoupage(x1)
    # On crée une liste vide
    x3=''
    sans_erreurs=[]
    avec_erreurs=[]
    # on parcourt la liste de liste "x2" obtenue par la fonction "decoupage"
    for i in x2:
        # le codage de i element de x2
        x4=codage(i, A)
        sans_erreurs.append(x4)
        # On convertit la liste obtenue après chiffre en entier pour pouvoir l'avoir une erreur
        x5=int(''.join(map(str, x4)), 2)
        # On génère la position de l'erreur allant de droite vers la gauche
        x6=ra.randint(1, 7)
        # On ajoute l'erreur
        x5=x5^(2**(x6 -1 ))
        avec_erreurs.append(map(int, np.binary_repr(x5, 7)))
        #avec_erreurs.append(map(int, list(np.binary_repr(x7, 7))))
        # On remplit x3 après avoir convertit le bloc chiffré en caractère en utilisant "chr"
        # NB: "chr" permet d'obtenir le caractère correspondant au code ascii donné
        x3+=chr(x5)
    return x1, x2, sans_erreurs, avec_erreurs, x3
#=========================================================================================
#
#
#                           PARTIE II: DECHIFFREMENT
#
#
#=========================================================================================
#
#=========================================================================================
#    On définit une fonction "decodage1" qui va decoder juste un seul caractère
# elle en entrée le caractère à décoder et la matrice de parité 
#=========================================================================================
def decodage1(a, H):
     # On récupère le code ascii de notre caractère
     x1=ord(a)
     # On le convertit en binaire sur 7 bits
     x2=np.binary_repr(x1, 7)
     # On transforme la chaine des bits en une liste d'entier constituée de ces bits
     x3=map(int, x2)
     # On crée un matrice de taille (1, 3) rempli de zero pour pouvoir evaluer le syndrome
     x4=np.zeros([1, 3], dtype=int)
     # On calcule le syndrome
     syndrome=np.dot(x3, np.transpose(H))
     syndrome2=[syndrome[i]%2 for i in xrange(3)]
     syndrome3=int(''.join(map(str, syndrome2)), 2)
     # On crée une liste vide x5 pour contenir le resultat
     x5=[]
     # On évalue le syndrome
     if syndrome3==0:
          x5.append(x3[:4])
     else:
          x3[syndrome3 -1]=x3[syndrome3 - 1]^1
          x5.extend(x3[:4])
     return x5
#=========================================================================================
# On définit une fonction "decodage2" qui utilise la fonction "decodage1" pour décoder tout
#  une cahine de caractère
#=========================================================================================
def decodage2(mot):
     x=mot
     H=matrice_de_control_de_hamm(3)
     x1=''
     x3=[]
     for i in x: 
          x2=decodage1(i, H)
          x3.extend(x2)
     x3=''.join(map(str, x3))
     x6=len(x3)
     for i in xrange(0, x6, 8):
         x1+=chr(int( x3[i:i+8], 2))
     return x1
#=========================================================================================
#
#
#
#
#
#
#=========================================================================================
#
#                     Fonction principale pour le chiffrement
#
#=========================================================================================
def fonction_principale():
    mot=raw_input("Saisir le mot à chiffrer:    \n")
    print('Veuillez suivre les indication suivantes et les utilisées dans l''ordre: \n')
    print('Tapez 1 pour convertir votre message en binaire \n')
    print('Tapez 2 pour découper votre message converti en binaire \n')
    print('Tapez 3 pour chiffrer votre les blocs \n')
    print('Tapez 4 pour ajouter des erreurs aux blocs chiffrés \n')
    print('Tapez 5 pour afficher le message chiffre français: \n')
    G=matrice_de_Gene_de_hamm(3)
    (m1, m2, m3, m4, m5)=codage2(mot, G)
    for i in xrange(5):
        x1=input()
        while x1 not in xrange(1, 6):
            print('Veuillez saisir une valeur comprise entre 1 et 4 \n')
            print('ou saisir drectement 1 pour voir le message chiffrer \n')
            x1=input()
        if x1==1:
            print('Votre message converti en binaire est: \n')
            print(m1)
        if x1==2:
            print('Votre découpé du message en bloc devient: \n')
            print(m2)
        if x1==3:
            print('Le chiffrement des blocs sans l''ajout des erreurs est \n')
            print(m3)
        if x1==4:
            print('Après l''ajout des erreurs, les blocs deviennent: \n')
            print(m4)
        if x1==5:
            print('Votre message chiffrée est: \n')
            print(m5)
    print('Voulez-vous déchiffrer le message ?')
    print('Veuillez saisir (O ou N): \n')
    reponse=raw_input()
    while reponse<>'O' and reponse<>'N':
        print('Voulez-vous déchiffrer le message ?')
        print('Veuillez saisir (O ou N): \n')
        reponse=raw_input()
    if reponse=='O':
        texteclair=decodage2(m5)
        print('Le message que vous avez chiffré est: \n')
        print(texteclair)
    print('Au revoir. Merci.............!!!!!!!!!!!!! \n')
