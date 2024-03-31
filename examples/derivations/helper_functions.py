import sympy as sp


####################################################
# Get i such that Ak[i] (or Dk[i]) appear in expr, #
# but Ak[i+1] (or Dk[i+1]) do not.                 #
####################################################
# (In our calculations we assume that once e.g. Ak[i+1] is not in expr 
#  anymore, neither is A[k+j] for any j = 1, 2, 3, ...; this holds true for
#  the calculations we are interested in.)

def get_maximal_D(expr):
   #
   str_expr = str(expr)
   #
   i = 1
   while ('D_{0}'.format(i) in str_expr):
      i += 1
   #
   return i - 1

def get_maximal_A(expr):
   #
   str_expr = str(expr)
   #
   i = 0
   while ('A_{0}'.format(i) in str_expr):
      i += 1
   #
   return i - 1


##################################################
# Helper functions for turning sympy expressions #
# into strings for saving to text files          #
##################################################

def single_sympy_expression_to_string(expr,
                                       N=100):
   #
   current_string = str(expr)
   #
   for i in range(0,N+1):
      i = N-i
      current_string = current_string.replace('A_{0}'.format(i),
                                 'Ak[{0}]'.format(i))
      current_string = current_string.replace('D_{0}'.format(i),
                                 'Dk[{0}]'.format(i))
   #
   current_string = current_string.replace('pi','sympy.pi')
   current_string = current_string.replace('exp','sympy.exp')
   current_string = current_string.replace('ln','sympy.ln')
   current_string = current_string.replace('log','sympy.ln')
   current_string = current_string.replace('\epsilon','epsilon')
   current_string = current_string.replace('x','xDL')
   current_string = current_string.replace('exDLp','exp')
   current_string = current_string.replace('sqrt','sympy.sqrt')
   #
   return current_string


def sympy_expression_to_string(expr):
    #
    # check if expression is a sum or a product
    if isinstance( expr, sp.core.mul.Add):
        list_of_terms = list(expr.args)
    else:
        list_of_terms = [expr]
    #
    list_of_strings = []
    for i,term in enumerate(list_of_terms):
        #
        string = single_sympy_expression_to_string(term)
        #
        if i > 0:
            #
            if string[0] != '-':
                string = '+' + string 
            #
        if i < len(list_of_terms) - 1:
            string = string + '\\\n   '
        list_of_strings.append(string)
    #
    string_out = ''
    for i,string in enumerate(list_of_strings):
        string_out = string_out + string 
    return string_out