import numpy as np
from scipy.linalg import expm

A = np.zeros((6, 6), dtype=complex)
for i in range(6):
    A[i, i] = 1.8 - 0.1 * i

ii = 1.2
for i in range(6):
    for j in range(i + 1, 6):
        A[i, j] = complex(ii, ii - 0.1)
        A[j, i] = complex(ii, -(ii - 0.1))
        ii -= 0.2
#print(np.trace(A))
A -= np.trace(A)*np.eye(6) / 6.0



A[0,0] = 0.0582178
A[1,1] = -0.0675741
A[2,2] = 0.0161615
A[3,3] = -0.0582178
A[4,4] = 0.0675741
A[5,5] = -0.0161615
A[0,1] = complex ( -0.0438132 , 0.0564889 )
A[0,2] = complex ( 0.0317717 , -0.0186406 )
A[0,3] = complex ( 0.0336684 , -0.028291 )
A[0,4] = complex ( -0.0264951 , 0.0536321 )
A[0,5] = complex ( -0.0776531 , 0.0864606 )
A[1,2] = complex ( 0.0650558 , -0.0119279 )
A[1,3] = complex ( -0.0126086 , 0.00982487 )
A[1,4] = complex ( -0.0586652 , -0.0984199 )
A[1,5] = complex ( 0.0306492 , -0.113321 )
A[2,3] = complex ( 0.0438132 , -0.0564889 )
A[2,4] = complex ( -0.0989423 , 0.00774407 )
A[2,5] = complex ( 0.0426485 , 0.0128435 )
A[3,4] = complex ( 0.0609364 , 0.00241003 )
A[3,5] = complex ( -0.0317717 , 0.0186406 )
A[4,5] = complex ( -0.0336684 , 0.028291 )

for i in range(6):
    for j in range(i + 1, 6):
        A[j,i] = complex(A[i,j].real, -A[i,j].imag)


#print(A)

def get_c(N, sign):
    c = np.zeros((2, N + 1), dtype=float)
    c[0, 0] = 1
    c[1, 0] = 1

    for k in range(1, N+1):
        c[0, k]=c[0, k-1] / k
        if k % 2 == 1:
            c[1, k] = -c[0, k]
        else:
            c[1, k]=c[0, k]
    return c[sign, :]

def get_p(A):
    A2 = np.dot(A, A)
    A3 = np.dot(A, A2)

    tr2 = np.trace(A2).real
    tr3 = np.trace(A3).real
    tr4 = np.trace(np.dot(A2, A2)).real
    tr5 = np.trace(np.dot(A2, A3)).real
    tr6 = np.trace(np.dot(A3, A3)).real

    print('Tr = ', [tr2, tr3, tr4, tr5, tr6])

    p = np.zeros(5)
    p[0] = -1.0/6.0*tr6 + 1.0/18.0*tr3*tr3 - 1.0/48.0*tr2*tr2*tr2 + 1.0/8.0*tr4*tr2
    p[1] = -1.0/5.0*tr5 + 1.0/6.0*tr2*tr3 
    p[2] = -1.0/4.0*tr4 + 1.0/8.0*tr2*tr2 
    p[3] = -1.0/3.0*tr3
    p[4] = -1.0/2.0*tr2
    
    return p

def get_q(p, c, N):
    q = np.zeros(6)
    if N > 5:
        ic = N - 6
        for i in range(6):
            q[i]=c[ic+1+i]
        # print('q = ', q)
        while ic >= 0:
            q5=q[5]
            q[5] = q[4]
            q[4] = q[3] -q5*p[4]
            q[3] = q[2] -q5*p[3]
            q[2] = q[1] -q5*p[2]
            q[1] = q[0] -q5*p[1]
            q[0] = c[ic]-q5*p[0]
            ic -= 1
            # print('q = ', q)
    else:
        for k in range(6):
            if k <= N:
                q[k] = c[k]
            else:
                q[k]=0.0
    return q

N = 30
c = get_c(N, 0)
print('c = ', c)

p = get_p(A) 
print('p = ', p)

q = get_q(p, c, N)
print('q = ', q)

A2 = np.linalg.matrix_power(A, 2)
A3 = np.linalg.matrix_power(A, 3)
A4 = np.linalg.matrix_power(A, 4)
A5 = np.linalg.matrix_power(A, 5)


expA = q[0]*np.eye(6) + q[1]*A + q[2]*A2 + q[3]*A3 + q[4]*np.dot(A2, A2) + q[5]*np.dot(A2, A3) 
for i in range(6):
    for j in range(6):
        print("{:10.8g} + {:10.8g}i \t".format(expA[i, j].real, expA[i, j].imag ), end='')
    print()

print()
print()

exact = expm(A)
for i in range(6):
    for j in range(6):
        print("{:10.8g} + {:10.8g}i \t".format(exact[i, j].real, exact[i, j].imag ), end='')
    print()

print(q)
