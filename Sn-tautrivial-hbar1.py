def degreelist(n,k):
	v = []
	if n == 1:
		v.append([k])
	else:
		r = []
		for a in range(0,k+1):
			r = degreelist(n-1,a)
			m = len(r)
			for b in range(0,m):
				r[b].append(k-a)
				v.append(r[b])
	return v

def degreematrix(n,k,p,b,oldDict,newDict,echelonmatrix):
        R = FractionField(PolynomialRing(GF(p), 'c'))
	f = []
	w = []
        d = degreelist(n,k)
	for q in b:
	    for y in b:
		f.append(beta(n,q,y,oldDict,newDict,R,echelonmatrix))
	    w.append(f)
	    f = []
	return Matrix(w)

def hilbertpolynomial(n,p): #S_n, char. p, l is a partition
    R = FractionField(PolynomialRing(GF(p),'c'))
    P = PolynomialRing(R,n,'x_')
    X = P.gens()
    generators = []
    oldDict = {}
    k = 0
    m = 1
    V.<t> = PolynomialRing(QQ)
    A = matrix([])
    s = 0
    echelonmatrixlist = [matrix([[]]),matrix([[]])]
    basislist = [[],[]]
    b = degreelist(n,k)
    while true:
	newDict = {}
        A = degreematrix(n,k,p,b,oldDict,newDict,echelonmatrixlist[0])
        m = rank(A)
        s = s + m*t^k
        A = basis(kernel(A)) #this is in truncated basis
        dOld = degreelist(n,k)
        for i in A:
            gen = 0
            for j in range(0,len(i)):
                if i[j] != 0:
                    monom = P(1)
                    for l in range(0,n):
                        monom = monom*(X[l]^(b[j][l]))
                    gen = gen + i[j]*monom
            print gen
        if m == 0:
            print s
            print generators
            return s, generators
	print m, "*t^",k
        A2 = [] 
        for x in A:
            a = []
            for i in range(0,len(dOld)):
                a.append(0)
            for y in range(0,len(x)):
                a[dOld.index(b[y])] = x[y]
            A2.append(a)
        A = A2
        if A != []:
            if echelonmatrixlist[1] != matrix([[]]):
                A = (matrix(A).stack(echelonmatrixlist[1])).rows()
        else:
            if echelonmatrixlist[1] != matrix([[]]):
                A = echelonmatrixlist[1].rows()
        d = degreelist(n,k+1)
        if len(A) != 0:
            B = []
            a = []
            for x in A:
               for i in range(0,n):
                   for e in range(0,len(d)):
                       a.append(0)
                   for y in dOld:
                       j = dOld.index(y)
                       z = list(y)
                       z[i] = z[i] + 1
                       w = d.index(z) 
                       a[w] = x[j]
                   B.append(a)
                   a = []
            B = matrix(basis(matrix(B).row_space())).echelon_form()
            echelonmatrixlist[0] = echelonmatrixlist[1]
            echelonmatrixlist[1] = B
            f = B.nrows()
            b = []
            toRemove = [] #get rid of all pivots
            for g in range(0,f):
                toRemove.append(first(B[g]))
            for i in range(0,len(d)):
                if not(i in toRemove):
                    b.append(d[i])
        else:
            b = d
        basislist[0] = basislist[1]
        basislist[1] = b   
        k = k + 1
	oldDict = newDict

def first(L_1): #spits out (one less than, since list indices start from 0) index of first non-zero entry of list
    for d in range(0,len(L_1)):
        if L_1[d]!=0:
            return d

def getAns(n,L_1,L_2,echelonmatrix,oldDict,R):
    if (n,tuple(L_1),tuple(L_2)) in oldDict:
	return oldDict[(n,tuple(L_1),tuple(L_2))]
    elif (n,tuple(L_2), tuple(L_1)) in oldDict:
        return oldDict[(n,tuple(L_2),tuple(L_1))]
    else:
        k = sum(L_1)
        d = degreelist(n,k)
        i_1 = d.index(L_1)
        j = echelonmatrix.ncols()
        m = echelonmatrix.nrows()
        piv = 0 # find pivot
        flag = true
        while flag:
            if piv > m-1:
                flag = false
            elif echelonmatrix[piv][i_1] != 0:
                flag2 = true
                for s in range(piv,i_1):
                    if echelonmatrix[piv][s] != 0:
                        flag2 = false
                if flag2:
                    flag = false
                else: 
                    piv = piv+1
            else:
                piv = piv+1
        if piv>m-1: #L_1 is not a pivot, so L_2 must be a pivot
            temp=getAns(n,L_2,L_1,echelonmatrix,oldDict,R)
            oldDict[(n,tuple(L_1),tuple(L_2))] = temp
            return temp
        f = R(0)
        for s in range(i_1 + 1, j):
            f = f - R(echelonmatrix[piv][s]*getAns(n,d[s],L_2,echelonmatrix,oldDict,R))
        oldDict[(n,tuple(L_1),tuple(L_2))] = R(f)
        return R(f)

def beta(n,L_1,L_2,oldDict, newDict, R,echelonmatrix): 
    c = R.gen(0)
    #n for S^n, list of exponents for x, then y, coupled with tau and tau* vectors, then partition p. to be consistent, sums of L_1 and L_2 are n, p is a partition of n, and v_1 and v_2 have the same dimension as representation for p.

#example: trivial representation over S_3, want beta(x1,y2). then n=3, first list is [1,0,0] for x1, and [0,1,0] for y2. the space is 1-dimensional so both v_1 and v_2 are (1), and then p=[3]. input is beta(3,[1,0,0],vector([1]),[0,1,0],vector([1]),[3]). note that before entering the inputs we have to check the dimension of S^(p) so that we have the right number of coordinates for v_1 and v_2.

#syntax - n is an integer, L_1 and L_2 are lists [], v_1 and v_2 are vectors ().
#in comments, ordinal numbers count up from 1, but in code, they count up from 0.
#    if (n,tuple(L_1),tuple(v_1),tuple(L_2),tuple(v_2),tuple(p)) in oldDict:
#        return oldDict[(n,tuple(L_1),tuple(v_1),tuple(L_2),tuple(v_2),tuple(p))]
    if sum(L_1)!=sum(L_2): #degrees not equal means beta=0
        return 0
    if sum(L_1)==sum(L_2)==0: #if the H half of the H\otimes\tau arguments are 1, we take the valuation tau*(tau), which is just the dot product when we use dual bases
        b = R(1)
        newDict[(n,tuple(L_1),tuple(L_2))] = b 
        return b
    else:
        b=R(0) #this will be the value of beta
        k=first(L_1) #finds the variable to throw over to the other side
        newL1 = list(L_1)
        newL1[k]=L_1[k]-1 #takes degree of (k+1)-th variable in 1st argument down by 1 (thrown over to the other side)
        for i in range(1,n):
            for j in range(i+1,n+1): #gives us all s_ij with i<j
                if k==i-1 and L_2[i-1]>L_2[j-1]: #(we have to shift indices) this deals with <y_k,alpha>: in this case it is 1
                    for q in range(0,L_2[i-1]-L_2[j-1]):
                        N_q=list(L_2)
                        N_q[i-1]=L_2[i-1]-1-q
                        N_q[j-1]=L_2[j-1]+q #each N_q corresponds to a summand in (prod(x)-prod(s(x)))/(alpha) term in the dunkl operator.
			b = b - c*getAns(n,newL1,N_q,echelonmatrix,oldDict,R)
                elif k==j-1 and L_2[i-1]>L_2[j-1]: #(we have to shift indices) this deals with <y_k,alpha>: in this case it is -1
                    for q in range(0,L_2[i-1]-L_2[j-1]):
                        N_q=list(L_2)
                        N_q[i-1]=L_2[i-1]-1-q
                        N_q[j-1]=L_2[j-1]+q #each N_q corresponds to a summand in (prod(x)-prod(s(x)))/(alpha) term in the dunkl operator.
			b = b + c*getAns(n,newL1,N_q,echelonmatrix,oldDict,R)
                elif k==i-1 and L_2[i-1]<L_2[j-1]: #everything is just reversed here
                    for q in range(0,L_2[j-1]-L_2[i-1]):
                        N_q=list(L_2)
                        N_q[i-1]=L_2[i-1]+q
                        N_q[j-1]=L_2[j-1]-1-q
			b = b + c*getAns(n,newL1,N_q,echelonmatrix,oldDict,R)
                elif k==j-1 and L_2[i-1]<L_2[j-1]: #everything is just reversed here
                    for q in range(0,L_2[j-1]-L_2[i-1]):
                        N_q=list(L_2)
                        N_q[i-1]=L_2[i-1]+q
                        N_q[j-1]=L_2[j-1]-1-q
			b = b - c*getAns(n,newL1,N_q,echelonmatrix,oldDict,R)
                else: #the sum cancels in this case
                    b=b
        if L_2[k] > 0:
            newL2=list(L_2)
            newL2[k]=newL2[k]-1 #takes the degree of the (k+1)-th variable in 2nd argument down by 1 to take the partial derivative
            b = b + (L_2[k])*getAns(n,newL1,newL2,echelonmatrix,oldDict,R)
        newDict[(n,tuple(L_1),tuple(L_2))] = b
        return b
