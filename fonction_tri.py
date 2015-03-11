def fonction_tri(l):
    l1=[]
    for i in l:
        if type(i)==list:
            i=fonction_tri(i)
            for k in i:
                l1.append(k)
        else:
            l1.append(i)
    return l1
            
