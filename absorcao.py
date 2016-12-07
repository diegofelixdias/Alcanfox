from sympy import *
a = Matrix(3, 3, [1, 4, 9, 7, 9, 3, 2, -2, 8])
print(a.is_diagonalizable())

(P, D) = a.diagonalize()
print(P)
