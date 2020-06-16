from sympy import Symbol
import numpy as np

# From the exercise:
def function_exercise(wxyz):
    w, x, y, z = wxyz
    return [(w**2)+x-(3*y)+(4*z)+(3/4),
            (3*w**2)+x-(y**2)+(z**2)+(13/4),
            (5*w)+(3*x**2)+y-(4*z**2)-(99/2),(8*w**2)-(14*x)+(6*y**2)-(7*z**2)+7]

def jacobian_exercise(wxyz):
    w, x, y, z = wxyz
    return [[2*w,1,-3,4],
            [6*w,1,-2*y,2*z],
            [5,6*x,1,-8*z],[16*w,-14,12*y,-14*z]]

def JF(x):
    
    J_save = np.array(jacobian_exercise(x))
    
    return J_save

def F(x):
    
    F_save =  np.array(function_exercise(x))
    
    return F_save

def iterative_newton(fun, var1, var2, jacobian,epsilon,itera):

    x_last = var1
    
    x_next = var2

    for i in range(itera):

        J = np.array(jacobian(x_last))
        
        F = np.array(fun(x_next))

        diff = np.linalg.solve( J, -F )
        
        x_last = x_last + diff
        
        x_next = x_next + diff

        # Validación del margen de error
        if np.linalg.norm(diff) < epsilon:
            
            return x_last
            
            break

    return x_last

# For the exercice:
x_sol = iterative_newton(function_exercise, [2.0,1.0,2.0,1.0],[2.0,1.0,2.0,1.0], jacobian_exercise, 0.1,100)
print('La solucion es:', x_sol )
       

def solve(x_0, tol, itera):

            a = (1) 

            b = (-2)

            for i in range(itera):

               Jx_kF =  iterative_newton(function_exercise, x_0, x_0,jacobian_exercise,tol,itera)

               y_k = x_0 - 1/2*(Jx_kF)

               Jy_kF = iterative_newton(function_exercise, y_k, x_0, jacobian_exercise,tol,itera)
            
               z_k = x_0 - Jy_kF

               Jy_xF = iterative_newton(function_exercise, y_k, x_0, jacobian_exercise,tol,itera)

               Jx_zF = iterative_newton(function_exercise, x_0, z_k, jacobian_exercise,tol,itera)

               Jy_zK = iterative_newton(function_exercise, y_k, z_k, jacobian_exercise,tol,itera)

               x_k = z_k - (a * Jx_zF) - (b * Jy_zK)

               x_0 = x_k
            
            return x_k
        
    
    
    

#b = [-3/4,-13/4,99/2,-7]
#x = 0
#A = [[1,1,-3,4],[3,1,-1,1],[5,3,1,-4],[8,-14,6,-7]]
    
            
    
    
