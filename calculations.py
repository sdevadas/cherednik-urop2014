from Sntautrivialhbar1 import hilbertpolynomial

calculations=[]

for p in [2,3,5,7,11]:
	n=p
	while n < 15:
		calculations.append((n,p,hilbertpolynomial(n,p)))
		n+=p
		
for l in calculations:
	print l
		