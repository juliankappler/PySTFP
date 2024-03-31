import sympy as sp
from math import factorial 


epsilon = sp.symbols(r'\epsilon',real=True,positive=True)

 # basic units of length, time, and diffusivity
L = sp.symbols('L',real=True,positive=True)
T = sp.symbols('T',real=True,positive=True)
D_0 = L**2 / T
D_x0 = sp.symbols(r'D(x_{0})',real=True,positive=True)

# other symbolic variables
dt = sp.symbols('\Delta{t}',real=True,positive=True)
dx = sp.symbols('\Delta{x}',real=True)
x = sp.symbols('x',real=True)
xDL = sp.symbols(r'\tilde{x}',real=True)
tDL = sp.symbols(r'\tilde{t}',real=True,nonnegative=True)
#

tD = sp.symbols(r'\tau_D',real=True,positive=True)

N_max = 10
Dk = sp.symbols('\mathcal{D}_0:%d'%N_max) # Eq. (29) of Ref. [1]
Ak = sp.symbols('\mathcal{A}_0:%d'%N_max) # Eq. (28) of Ref. [1]

dA = sp.symbols(r'a_0:%d'%N_max) # derivatives of a, evaluated at x0
dD = [sp.symbols(r'd_0',real=True,positive=True)]
for i in range(1,N_max):
    dD.append(sp.symbols(r'd_{0}'.format(i),real=True))

R = sp.sqrt( 2*dD[0]*dt )
epsilon = R / L
tauD = L**2 / dD[0]




# When we retrieve symbolic expressions from a PySTFP instance, those 
# expressions use symbols that are defined within the instance.
# We substitute these symbols with the symbols defined in the cells above
def substitute_local_variables_in_expr(
                expr, # sympy expression
                remote, # namespace of the symbols that we want to replace
                    ):
    #
    expr = expr.subs(remote.L,L)
    expr = expr.subs(remote.T,T)
    expr = expr.subs(remote.t,dt)
    expr = expr.subs(remote.x,x)
    expr = expr.subs(remote.xDL,xDL)
    expr = expr.subs(remote.tDL,tDL)
    expr = expr.subs(remote.epsilon,epsilon)
    for k in range(N_max):
        expr = expr.subs(remote.Dk[k],Dk[k])
        expr = expr.subs(remote.Ak[k],Ak[k])
    return expr 

def substitute_local_variables_in_dictionary(
                dictionary,
                remote,
            ):
    #
    output_dictionary = {}
    #
    for i, expr in dictionary.items():
        if isinstance(expr, sp.Expr):
            output_dictionary[i] = substitute_local_variables_in_expr(
                                            expr=expr,
                                            remote=remote,
                                            )
        else:
            output_dictionary[i] = expr
    #
    return output_dictionary






# Functions for substituting the definitions of \mathcal{A}, \mathcal{D},
# c.f. Eqs. (28), (29) of Ref. [1]

def get_Dk(n):
    if n == 0:
        return 1
    D_n = L**n/dD[0]/factorial(n) * dD[n]
    return D_n

def get_Ak(n):
    A_n = L**(n+1)/dD[0]/factorial(n) * dA[n]
    return A_n


def substitute_physical_dimensions_in_expr(expr):
    #
    for k in range(N_max):
        expr = expr.subs(Dk[k],get_Dk(k))
        expr = expr.subs(Ak[k],get_Ak(k))
    #
    expr = expr.subs(xDL,dx/R)
    expr = expr.subs(tDL,dt/tauD)
    #
    return expr

def substitute_physical_dimensions_in_dictionary(dictionary):
    #
    dictionary_output = {}
    #
    for i, expr in dictionary.items():
        if isinstance(expr, sp.Expr):
            dictionary_output[i] = substitute_physical_dimensions_in_expr(
                                            expr=expr
                                            )
        else:
            dictionary_output[i] = expr
    #
    return dictionary_output




def substitute_physical_dimensions(
                        remote_dictionary,
                        remote,
                            ):
    #
    dictionary_local = substitute_local_variables_in_dictionary(
                dictionary=remote_dictionary,
                remote=remote,
                    )
    #
    dictionary_output = substitute_physical_dimensions_in_dictionary(
                    dictionary=dictionary_local
                        )
    #
    return dictionary_output
