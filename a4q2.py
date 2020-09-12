#extract matrix info
with open('a4q2mat.txt','r') as f:
	a =[[float(num) for num in line.split(',')] for line in f]


n = len(a)
#set l and u to zero matrices of nxn
l = [[float(0) for x in range(n)]
			for y in range(n)]
u =[[float(0) for x in range(n)]
			for y in range(n)]
I =[[float(0) for x in range(n)]
			for y in range(n)]
I2 =[[float(0) for x in range(n)]
			for y in range(n)]
I3 =[[float(0) for x in range(n)]
			for y in range(n)]


#function fforr partial pivoting
def partialpivoting(a):
	for i in range(0,n-1):
		if abs(a[i][i]) == 0 :
			for j in range(i,n):
				if abs(a[j][i])>abs(a[i][i]) :
					a[i],a[j]=a[j],a[i]
					
		else:
			pass			
	return a

a=  partialpivoting(a)

#function to get cofactors
def cofactor(a,t,p,q,n):
	i=0
	j=0

	for k in range(n):
		for l in range(n):
			if k != p and l != q:
				t[i][j] = a[k][l]
				j += 1

				if j== n-1:
					j=0
					i+=1

#function to get determinant using cofactor
def determinant(a,n):
	D =0
	if n ==1:
		return a[0][0]
	#define a temp matrix to form cofactor
	t =[[float(0) for x in range(n)]
			for y in range(n)]	
	sign = 1

	for f in range(n):
		cofactor(a,t,0,f,n)
		D += (sign*a[0][f]*determinant(t,n-1))
		#to alternate sign
		sign = -sign
	return D
#function to check det(a)  
def invertable(a,n):
	if determinant(a,n) != 0:
		return True
	else:
		return False

#output
if invertable(a,n):
	print(f'yes, a is invertable')
else:
	print(f'non invertable')

#using Doolittle algorithm
def LUdecomposition(a,l,u):

	for i in range(n):
		for j in range(n):
				sum1 = 0
				for k in range(n):
					sum1 += (l[i][k]*u[k][j])
				#using conditionals to get L and U 
				if i == j :
					l[i][j] = 1

					u[i][j]= a[i][j] - sum1
				elif i <j:
					u[i][j]= a[i][j] - sum1
				else :
				
					l[i][j]=float((a[i][j]-sum1)/u[j][j])	
					
			

			
	return a,l,u
A,L,U = LUdecomposition(a,l,u)



#define gaussian jordan functtion 
def gaussianjordan(a,I):
	#define pivot and loop over the firrrst list
	
	for i in range(0,n):

		pivot = a[i][i]
		#make diagonals of zero matrix 1 to get identity matrix
		I[i][i]=1

		#loop overr second list to trravvel forrr rrow operation
		for j in range(0,n):
			a[i][j] = a[i][j]/pivot
			I[i][j] = I[i][j]/pivot

		for k in range(0,n):
			factor = a[k][i]
			if k == i or factor== 0:
				pass
			
			else:
				
				for j in range(0,n):
					a[k][j] =a[k][j]-factor*a[i][j]
					I[k][j] =I[k][j]-factor*I[i][j]

			#else:
			#continue
					
	return I
#tthus we gett invverrse of L and U
X= gaussianjordan(L,I)



Y= gaussianjordan(U,I2)


#usingg A^-1 = U^-1 x L^-1 
def crossproduct(a,b)	:
	cross =[[float(0) for x in range(n)]
			for y in range(n)]
	
	for i in range(0,n):
		for j in range(0,n):
			for k in range(0,n):
				
				cross[i][j] += a[i][k]*b[k][j]

	return cross
C = crossproduct(Y,X)
print(f'Inverse of A (using LU): {C}')


Z = gaussianjordan(A,I3)
print(f'Inverse of A (using gaussian jordan):{Z}')


'''
yes a is invertable
Inverse of A (using LU): [[0.3333333333333333, -0.2500000000000001, 1.666666666666667, -1.8333333333333333], [0.0, 0.08333333333333337, -0.666666666666667, 0.8333333333333333], [0.0, 0.16666666666666666, -0.33333333333333326, -0.3333333333333333], [0.0, -0.08333333333333333, 0.6666666666666666, 0.16666666666666666]]
Inverse of A (using gaussian jordan):[[0.3333333333333333, -0.2500000000000001, 1.666666666666667, -1.8333333333333333], [0.0, 0.08333333333333337, -0.666666666666667, 0.8333333333333333], [0.0, 0.16666666666666666, -0.33333333333333326, -0.3333333333333333], [0.0, -0.08333333333333333, 0.6666666666666666, 0.16666666666666666]]

'''
