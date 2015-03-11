def nbr_premier():
    for m in range(2,10):
        for i in range(2,m):
            if m % i == 0:
                print m,"n'est pas premier" , "m = ",i," x " ,m/i
                break
        else: # a executer si on sort de la boucle sans un break
            print m, "est premier"
