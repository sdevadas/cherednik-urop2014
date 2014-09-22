from Sntautrivialhbar1 import hilbertpolynomial

calculations=[]

nplist = [(2,2),(4,2),(3,3),(6,3),(6,2),(9,3), (5,5)]

for n,p in nplist:
		calculations.append((n,p,hilbertpolynomial(n,p)))
		
for l in calculations:
	print l
		