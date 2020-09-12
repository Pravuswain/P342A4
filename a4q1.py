#extract matrix info
with open('a4q1mat1.txt','r') as f:
	a =[[float(num) for num in line.split(',')] for line in f]

with open('a4q1mat2.txt','r') as f:
	b =[float(num) for num  in f]

n = len(b)
#set l and u to zero matrices of nxn
l = [[float(0) for x in range(n)]
			for y in range(n)]
u =[[float(0) for x in range(n)]
			for y in range(n)]
Y =[float(0) for y in range(n)]
X =[float(0) for y in range(n)]
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

#solving equations using L and U
def ForwardBackwardSubstitution(L,U,b,Y,X):

#forward sustitution on Ly=b
	for i in range(0,n):
		sum2 = 0
		for j in range(0,i+1):
			
			sum2 += (L[i][j]*Y[j])
		Y[i]=b[i]-sum2
#backward substitution on Ux=y	
	for i in range(0,n):
		sum3 = 0
		for j in range(0,i+1):
			
			sum3 += (U[i][j]*X[j])
		X[i]=float((Y[i]-sum3)/U[i][i])

	return Y,X
y,x = ForwardBackwardSubstitution(L,U,b,Y,X)

print(f'LU decomposition of A gives:')
print(f"L={L}")
print(f'U={U}')

print(f"solutions of linear equations : ")			
for i in range(n):

	print(f' x{i+1}={x[i]}')


'''OUPUT :

LU decomposition of A gives:
L=[[1, 0.0, 0.0, 0.0], [0.0, 1, 0.0, 0.0], [1.0, 2.0, 1, 0.0], [2.0, 1.0, 1.5, 1]]
U=[[1.0, 0.0, 1.0, 2.0], [0.0, 1.0, -2.0, 0.0], [0.0, 0.0, 2.0, -2.0], [0.0, 0.0, 0.0, -2.0]]
solutions of linear equations : 
 x1=6.0
 x2=-3.0
 x3=-1.0
 x4=3.0


'''