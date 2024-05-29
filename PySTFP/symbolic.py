#!/usr/bin/env python

import numpy as np
import sympy
import time

import PySTFP.sympy_definitions.normalization_preserving_propagator as npp
import PySTFP.sympy_definitions.positivity_preserving_propagator as ppp
import PySTFP.sympy_definitions.moments as moments
import PySTFP.sympy_definitions.dS_m_dt as dS_m_dt
import PySTFP.sympy_definitions.dS_tot_dt as dS_tot_dt
import PySTFP.sympy_definitions.S_Gibbs as S_Gibbs


'''
Bibliography

[1] Short-time Fokker-Planck propagator beyond the Gaussian approximation,
    Julian Kappler
	arXiv: http://arxiv.org/abs/2405.18381

''';



class PySTFP:
	'''
	This class contains perturbative results for the short-time Fokker-Planck
	propagator

	'''
	def __init__(self,parameters={}):
		'''Initialization method of the class

		Sets simulation parameters and defines symbolic length-, time-,
		and diffusivity scales

		Args:
			parameters (dict): A dictionary with the parameters defined by the
							   user. This is optional, all parameters can be
							   also set after initialization of the class.

		Returns:
			None
		'''
		#
		# basic units of length, time, and diffusivity
		self.L = sympy.symbols('L',real=True,positive=True)
		self.T = sympy.symbols('T',real=True,positive=True)
		self.D0 = self.L**2 / self.T
		#
		# other symbolic variables
		self.t = sympy.symbols('t',real=True,nonnegative=True)
		self.x = sympy.symbols('x',real=True)
		self.xDL = sympy.symbols(r'\tilde{x}',real=True)
		self.tDL = sympy.symbols(r'\tilde{t}',real=True,nonnegative=True)
		self.epsilon = sympy.symbols(r'\epsilon',real=True,positive=True)
		#
		self.N_max = 11 # Instead of using N_max = 11, it would be preferable 
		# to find the largest symbol Ak[i], Dk[i] we need to have to 
		# accomodate all the pre-calculated perturbative results 
		self.N_max = self.__get_N_max()
		self.Dk = sympy.symbols('D_0:%d'%self.N_max)
		self.Ak = sympy.symbols('A_0:%d'%self.N_max)
		#
		#
		self.verbose = True
		#
		self.a = None
		self.D = None
		self.x0 = None
		#
		self.A_functions_constructed = False
		self.D_functions_constructed = False
		#
		self.factorials = []
		for i in range(self.N_max):
		    self.factorials.append(np.math.factorial(i))
		#
		# get perturbation series and set local variables
		self.Q_dict = {}
		for i,term in (npp.Q_dict).items():
			self.Q_dict[i] = self.__substitute_local_variables(
							expression=term,
							remote=npp
							)
		#
		self.P0 = self.__substitute_local_variables(
							expression=npp.P0,
							remote=npp
							)
		#
		self.Q_dict_ppp = {} # ppp = positivity-preserving propagator
		for i,term in (ppp.Q_dict).items():
			self.Q_dict_ppp[i] = self.__substitute_local_variables(
							expression=term,
							remote=ppp
							)
		#
		# get moments series and set local variables
		self.moments = {}
		for i,power_series_terms in (moments.moments_dict).items():
			self.moments[i] = {}
			for j,term in (power_series_terms).items():
				self.moments[i][j] = self.__substitute_local_variables(
								expression=term,
								remote=moments
								)
		#
		# get Sm_dot (medium entropy production rate) series and 
		# set local variables
		self.dS_m_dt = {}
		for i,term in (dS_m_dt.dS_m_dt_dict).items():
			self.dS_m_dt[i] = self.__substitute_local_variables(
								expression=term,
								remote=dS_m_dt
								)

		# get S_tot_dot (total entropy production rate) series and 
		# set local variables
		self.dS_tot_dt = {}
		for i,term in (dS_tot_dt.dS_tot_dt_dict).items():
			self.dS_tot_dt[i] = self.__substitute_local_variables(
								expression=term,
								remote=dS_tot_dt
								)

		# get S_Gibbs (Gibbs entropy) series and set local variables
		self.S_Gibbs = {}
		for i,term in (S_Gibbs.S_Gibbs_dict).items():
			self.S_Gibbs[i] = self.__substitute_local_variables(
								expression=term,
								remote=S_Gibbs
								)



	def set_parameters(self,parameters):
		'''Sets parameters of an instance of the class

		With this method the parameters of an existing instance of the class
		can be set or changed.

		Args:
			parameters (dict): A dictionary with the parameters defined by the
							   user.

		Returns:
			None
		'''
		#
		try:
			self.set_a(a=parameters['a'])
		except KeyError:
			pass
		#
		try:
			self.set_D(D=parameters['D'])
		except KeyError:
			pass
		#




	def check_dimensions(self,
					quantity, # quantity one wants to check 
					dimensions, # units one wants to check for
					return_with_dimensions=True,
					):
		''' Check if a variable is a dimensionless, or has given units

		Args:
			quantity: The variable of which one wants to check the units
			dimensions: The units one would like to check for
			return_with_dimensions: If True, then the function returns the
			                        quantity as a sympy expression with the
									dimensions. This is the default
									behavior.
									If False, then the function returns the
									quantitiy as a float, without its 
									dimensions.
		
		Returns:
			the quantity, depending on the arguments of the function call 
			either with units (as a sympy expression) or without units (as a 
			float).

		''';
		#
		try: # check if the quantity is a float
			quantity = float(quantity) 
		except TypeError: # if it is not a float
			try: # check if the quantity has units "dimensions"
				quantity = float(sympy.simplify(quantity/dimensions))
			except TypeError: # if it does not have the right units
				raise RuntimeError(
					"Provided quantity = {0} is neither".format(quantity) \
					+ " dimensionless nor in units of {0}.".format(dimensions)
						)
		#
		if return_with_dimensions:
			return quantity * dimensions
		else:
			return quantity 




	def set_a(self,a):
		'''Set drift

		This method sets the drift profile for both analytical and numerical
		functionalities of the the class.

		Args:
			a (sympy expression): Drift profile as function of self.x, self.t,
								  and in units self.T/self.T.
								  Example:
								    a = self.T/self.T * \
								  		sympy.sin(self.x/self.L - self.t/self.T)

		Returns:
			None
		'''
		self.a = a # in units of self.L/self.T
		self.A_functions_constructed = False
		#
		expr = sympy.expand( a / (self.L/self.T) )
		expr = expr.subs(self.x,self.x*self.L)
		expr = sympy.expand( expr )
		self.a_lambda = sympy.lambdify(self.x,expr)

	def set_D(self,D):
		'''Set diffusivity

		This method sets the diffusivity profile for both analytical and
		numerical functionalities of the the class.

		To be consistent with the assumptions of the Fokker-Planck equation,
		the diffusivity should always be strictly positive.

		Args:
			D (sympy expression): Diffusivity profile as function of self.x,
								  self.t, and in units
								  self.D0 = self.L**2 / self.T
								  Example:
								    D = self.D0*(1+0.1*sympy.sin(self.x/self.L))

		Returns:
			None
		'''
		self.D = D
		#
		self.D_functions_constructed = False
		#
		expr = D/self.D0
		expr = expr.subs(self.x,self.x*self.L)
		expr = sympy.expand( expr )
		#
		self.D_lambda = np.vectorize( sympy.lambdify(self.x,expr) )


	def set_x0(self,
				x0, # initial position for propagator
				evaluate_DA=True,
				):
		''' Set initial condition for the propagator

		Args:
			x0 (sympy.expr or float): Initial condition for the propagator. 
				If given as sympy.expr, has to have unit self.L. 
				If given as float, is assumed to be in units of self.L
			evaluate_DA (bool, default True): If True, then drift and 
				diffusivity will be evaluated at x0 and the results will be
				stored in the class instance. 
		
		Returns:
			None

		''';
		#
		self.x0 = self.check_dimensions(quantity=x0,
									dimensions=self.L,
									return_with_dimensions=False)
		#
		if evaluate_DA:
			self.__eval_Ak(x0=self.x0)
			self.__eval_Dk(x0=self.x0)
			self.__eval_D_x0(x0=self.x0)
		else:
			self.D_evaluated=False
			self.A_evaluated=False
			self.Dx_0_evaluated=False

	def set_D_x0(self,D_x0):
		''' Set value of D(x0)
		
		''';
		self.D_x0_float = self.check_dimensions(quantity=D_x0,
									dimensions=self.D0,
									return_with_dimensions=False)
		self.D_x0 = self.Dx_0_float * self.D0

	def construct_A_functions(self):
		''' Construct lambda functions for A_k from symbolic drift profile

		In Ref. [1] the A_k are written \tilde{\mathcal{A}}_k and defined in 
		Eq. (28).

		This method defines lists:
		- self.Ak_symbolic: sympy expressions for A_k, and 
		- self.Ak_lambda: lambda functions A_k(x0)

		Args: None
		Returns: None
		'''
		#
		error_msg_template = ["Cannot construct {0} functions: ",
					 "No {0} profile has been set. ",
					 "Please set a {0} profile via the class ",
					 "method {0}."]
		#
		if self.a is None:
			error_msg = error_msg_template[0].format('A_k') \
						+ error_msg_template[1].format('drift') \
						+ error_msg_template[2].format('drift') \
						+ error_msg_template[3].format('set_a')
			raise RuntimeError(error_msg)
		#
		if self.D is None:
			error_msg = error_msg_template[0].format('D_k') \
						+ error_msg_template[1].format('diffusivity') \
						+ error_msg_template[2].format('diffusivity') \
						+ error_msg_template[3].format('set_D')
			raise RuntimeError(error_msg)
		#
		self.Ak_symbolic = []
		#
		for i in range(self.N_max):
			#
			if i == 0:
				current_derivative = self.a
			else:
				current_derivative = sympy.diff(current_derivative,
												self.x,1)
			#
			self.Ak_symbolic.append( \
							  self.L**(i+1) / self.D \
							* sympy.Rational(1,self.factorials[i]) \
							* current_derivative \
										)
		#
		self.Ak_lambda = []
		for i, current_function in enumerate(self.Ak_symbolic):
			self.Ak_lambda.append(
				sympy.lambdify( self.x,
					sympy.simplify(current_function.subs(self.x,
												self.x*self.L) ))
									)
		#
		self.A_functions_constructed = True


	def construct_D_functions(self):
		'''Construct lambda functions for D_k from symbolic diffusivity profile

		In Ref. [1] the D_k are written \tilde{\mathcal{D}}_k and defined in 
		Eq. (29).

		This method defines lists:
		- self.Dk_symbolic: sympy expressions for D_k, and 
		- self.Dk_lambda: lambda functions D_k(x0)

		Args: None
		Returns: None
		'''
		#
		if self.D is None:
			raise RuntimeError("Cannot construct D functions: " \
					+ "No diffusivity profile has been set. " \
					+ "Please set a diffusivity profile via the class "\
					+ "method set_D.")
		#
		self.Dk_symbolic = []
		#
		for i in range(self.N_max):
			#
			if i == 0:
				current_derivative = self.D
			else:
				current_derivative = sympy.diff(current_derivative,
												self.x,1)
			#
			self.Dk_symbolic.append( \
							  self.L**i / self.D \
							* sympy.Rational(1,self.factorials[i]) \
							* current_derivative \
										)
		#
		self.Dk_lambda = []
		for i, current_function in enumerate(self.Dk_symbolic):
			self.Dk_lambda.append(
				sympy.lambdify( self.x,
						sympy.simplify(current_function.subs(self.x,
												self.x*self.L) ))
									)
		#
		self.D_lambda = sympy.lambdify(self.x, 
					sympy.simplify(self.D.subs(self.x,self.x*self.L)/self.D0),
											)
		#
		self.D_functions_constructed = True


	def __pad_array(self,array):
		'''

		to do


		'''
		# The arrays Dk, Ek need to be of length self.N_max. If they are
		# longer, we trim and throw a warning. If they are shorter, we
		# pad with zeros and throw a warning
		if len(array) > self.N_max:
			array_ = array[:self.N_max]
		elif len(array) < self.N_max:
			array_ = np.zeros(self.N_max,dtype=float)
			array_[:len(array)] = array
		else:
			array_ = array
		return array_




	def get_Dk_Ak_from_derivatives(self,
					D_derivatives,
					a_derivatives=None):
		''' Get dimensionless coefficients A_k, D_k from an array of 
		derivatives of drift and diffusivity

		In Ref. [1] the A_k, D_k are written \tilde{\mathcal{A}}_k, 
		\tilde{\mathcal{D}}_k. This method evaluates the right-hand side of
		Eqs. (28), (29) to calculate A_k, D_k.

		Args:
			D_derivatives (np.array): spatial derivatives of D at x0, such 
				that D_derivatives[n] = (\partial_x^n D)(x_0)
			a_derivatives (np.array): spatial derivatives of a at x0, such 
				that a_derivatives[n] = (\partial_x^n a)(x_0)
		Returns:
			Dk (np.array): coefficients Dk[n] = \tilde{\mathcal{D}}_n
			Ak (np.array): coefficients Ak[n] = \tilde{\mathcal{A}}_n

		''';
		#
		if a_derivatives is None:
			a_derivatives = np.zeros_like(D_derivatives)
		#
		n = min([len(D_derivatives),len(a_derivatives)])
		if n == 0:
			raise RuntimeError('Please provide non-empty arrays for the'\
				+' derivatives')
		#
		L = 1. # since the arrays with the derivatives are in units of self.L
		tD = L**2 / D_derivatives[0] 
		#
		Dk = np.zeros(n,dtype=float)
		Ak = np.zeros(n,dtype=float)
		#
		for i in range(n):
			Dk[i] = L**i * D_derivatives[i]/D_derivatives[0]/self.factorials[i]
			Ak[i] = tD * L**(i-1) / self.factorials[i] * a_derivatives[i]
		#
		return Dk, Ak
			

	def substitute_Dk_Ak(self,
						expr,
						Dk,
						Ak,
						):
		''' Substitute coefficients D_k, A_k into a sympy expression

		In Ref. [1] the D_k, A_k are written \tilde{\mathcal{D}}_k, 
		\tilde{\mathcal{A}}_k, and defined in Eqs. (28), (29).

		Args:
			expr (sympy.expr): expression in which we want to substitute the
				coefficients
			Dk (np.array): coefficients Dk[n] = \tilde{\mathcal{D}}_n
			Ak (np.array): coefficients Ak[n] = \tilde{\mathcal{A}}_n
		Returns:
			expr (sympy.expr): input expression with coefficients substituted

		''';
		#
		# pad arrays
		Ak = self.__pad_array(Ak)
		Dk = self.__pad_array(Dk)
		#
		for k, D_k in enumerate(Dk):
			expr = expr.subs(self.Dk[k],D_k)
		for k, A_k in enumerate(Ak):
			expr = expr.subs(self.Ak[k],A_k)
		#
		return expr


	def substitute_physical_time(self,
					expr,
					D_x0=None,
					):
		''' Substitute physical time into a sympy expression (that contains
		the dimensionless time variable)

		In Ref. [1] the relation between dimensionless and physical time is
		given by Eq. (13).

		Args:
			expr (sympy.expr): expression in which we want to substitute the
				physical time
			D_x0 (sympy.expr, float, or None): diffusivity at inital point x_0
				If None (default): Uses internally stored value for D(x_0)
				If sympy.expr: D(x_0), with units L^2/T
				If float: D(x_0), expressed in units of L^2/T
		Returns:
			expr (sympy.expr): input expression, but with physical time instead
				of dimensionless time

		''';
		#
		if D_x0 is None:
			D_x0 = self.D_x0
		else:
			D_x0 = self.check_dimensions(quantity=D_x0,
								dimensions=self.D0)
		time_conversion_factor = float(sympy.simplify(D_x0/self.D0))
		# convert time from units of tau_D to units of T
		return expr.subs(self.tDL, self.tDL * time_conversion_factor)

	def substitute_physical_position(self,
					expr,
					x0=None,
					):
		''' Substitute physical position into a sympy expression (that contains
		the dimensionless position variable)

		In Ref. [1] the relation between dimensionless and physical time is
		given by Eq. (14).

		Args:
			expr (sympy.expr): expression in which we want to substitute the
				physical position
			x0 (sympy.expr, float, or None): inital point x_0
				If None (default): Uses internally stored value for x_0
				If sympy.expr: x_0, with units L
				If float: x_0, expressed in units of L
		Returns:
			expr (sympy.expr): input expression, but with physical position 
				instead of dimensionless position

		''';
		#
		if x0 is None:
			x0 = self.x0
		else:
			x0 = self.check_dimensions(quantity=x0,
								dimensions=self.L)
		x0 = float(sympy.simplify(x0/self.L))
		#
		epsilon = sympy.sqrt(2 * self.tDL)
		#
		return expr.subs(self.xDL,
			(self.xDL - x0)/epsilon \
				)

		
	#########################
	# Short-time propagator #
	#########################

	def get_symbolic_perturbative_propagator(self,
						order,
						positivity_preserving=False,
						):
		''' Return symbolic expression for perturbative short-time propagator

		Returns either Eq. (23) or Eq. (41) from Ref. [1], in dimensionless
		units

		Args:
			order (int): order up to which the propagator should be returned. 
				Inclusive, i.e. the term K = order is also included
			positivitiy_preserving (bool): 
				if False (default): Return normalization-preserving propagator
					Eq. (23)
				if True: Return positivity-preserving propagator Eq. (41)

		Returns:
			expr (sympy.expr): symbolic short-time propagator

		''';
		#
		if positivity_preserving:
			Q_dict = self.Q_dict_ppp
		else:
			Q_dict = self.Q_dict
		#
		if len(Q_dict) <= order:
			raise RuntimeError("Perturbative solution not available to"\
				+ " requested order = {0}.".format(order) \
				+ " Please use order <= {0} ".format(len(Q_dict)-1))
		#
		epsilon = sympy.sqrt(2*self.tDL)
		perturbation_series = 0
		for k in range(order+1):
			perturbation_series += epsilon**k * Q_dict[k]
		#
		if positivity_preserving:
			P = self.P0 * sympy.exp(perturbation_series)
		else:
			P = perturbation_series * self.P0
		#
		return P


	def get_probability_density_DA(self,
				Dk, # numpy arrays, dimensionless quantities
				Ak, # numpy arrays, dimensionless quantities
				x0=None, # float, or explicitly in units of self.L
				D_x0=None, # float, or explicitly in units of self.D0
				order=4,
				dimensionless_position=False,
				dimensionless_time=False,
				positivity_preserving=False,
				):
		''' Return perturbative short-time propagator as lambda function

		Returns either Eq. (23) or Eq. (41) from Ref. [1], either with physical
		dimensions or in dimensionless units

		Args:
			Dk (np.array): coefficients Dk[n] = \tilde{\mathcal{D}}_n
			Ak (np.array): coefficients Ak[n] = \tilde{\mathcal{A}}_n
			x0 (sympy.expr, float, or None): inital point x_0
				If None (default): Uses internally stored value for x_0
				If sympy.expr: x_0, with units L
				If float: x_0, expressed in units of L
			D_x0 (sympy.expr, float, or None): diffusivity at inital point x_0
				If None (default): Uses internally stored value for D(x_0)
				If sympy.expr: D(x_0), with units L^2/T
				If float: D(x_0), expressed in units of L^2/T
			order (int): order up to which the propagator should be returned. 
				Inclusive, i.e. the term K = order is also included
			dimensionless_position (bool): 
				if False (default): return propagator using in physical lengths
				if True: return propagator using dimensionless distances
			dimensionless_time (bool): 
				if False (default): return propagator using in physical time
				if True: return propagator using dimensionless time
			positivitiy_preserving (bool): 
				if False (default): Return normalization-preserving propagator
					Eq. (23)
				if True: Return positivity-preserving propagator Eq. (41)

		Returns:
			P (lambda function): short-time propagator as a function of (x,t).
				Depending on the input arguments, (x,t) can be dimensionless
				(as defined in Eqs. (13), (14) of Ref. [1]) or in units of 
				length L and/or time T.

		''';
		#
		if x0 is None:
			x0 = self.x0
		#
		# recall that tDL is time, expressed in units of tau_D
		epsilon = sympy.sqrt(2 * self.tDL)
		#
		# get perturbative solution
		P = self.get_symbolic_perturbative_propagator(order=order,
						positivity_preserving=positivity_preserving)
		#
		# substitute numerical values for Ak, Dk
		P = self.substitute_Dk_Ak(expr=P,Dk=Dk,Ak=Ak)
		#
		if not dimensionless_position:
			# convert to length in units of L
			P = P / epsilon 
			P = self.substitute_physical_position(expr=P,x0=x0)
			#
		if not dimensionless_time:
			# convert to time in units of T
			P = self.substitute_physical_time(expr=P,D_x0=D_x0)
			#
		#
		P = sympy.simplify(P)
		#
		return sympy.lambdify((self.xDL,self.tDL),P)

	def get_probability_density_from_derivatives(self,
			D_derivatives, # numpy array, 
			  # D_derivatives[i] in units of self.D0 / self.L**i
			a_derivatives, # numpy array,
			  # a_derivatives[i] in units of self.L**(1-i)/self.T
			x0=None, # float, in units of self.L
			order=4,
			dimensionless_position=False,
			dimensionless_time=False,
			positivity_preserving=False,
			):
		''' Return perturbative short-time propagator as lambda function

		This is a convenience function that calls
			self.get_probability_density_DA.
		
		The only difference is that the present function accepts the 
		derivatives of a(x) and D(x) at x_0 as input arguments. It then 
		calculates the corresponding dimensionless coefficients D_k, A_k
		and calls self.get_probability_density_DA with those.

		''';
		#
		if x0 is None:
			x0 = self.x0
		#
		Dk, Ak = self.get_Dk_Ak_from_derivatives(
					D_derivatives=D_derivatives,
					a_derivatives=a_derivatives
						)
		#
		return self.get_probability_density_DA(
				Dk=Dk,
				Ak=Ak,
				order=order,
				dimensionless_position=dimensionless_position,
				dimensionless_time=dimensionless_time,
				positivity_preserving=positivity_preserving,
				x0=x0,
				D_x0=D_derivatives[0],
				)


	def get_probability_density(self,
				x0=None,
				order=2,
				positivity_preserving=False,
				dimensionless_position=False,
				dimensionless_time=False,
				):
		''' Return perturbative short-time propagator as lambda function

		This is a convenience function that calls
			self.get_probability_density_DA.
		
		The only difference is that the present function automatically uses
		the dimensionless coefficients D_k, A_k that are stored in the class
		instance. Using those, it calls self.get_probability_density_DA.

		''';
		#
		if x0 is not None:
			self.set_x0(x0=x0,
						evaluate_DA=True)
		#
		return self.get_probability_density_DA(
				Dk=self.Dk_evaluated,
				Ak=self.Ak_evaluated,
				order=order,
				dimensionless_position=False,
				dimensionless_time=False,
				positivity_preserving=positivity_preserving,
				)


	'''
	In the following section the perturbative observables (moments, entropy)
	are defined. I should add comments to the functions!
	''';

	############################
	# Moments < (\Delta x)^n > #
	############################

	def get_moment(self,
				x0=None,
				n=0,
				order=4,
				dimensionless_position=False,
				dimensionless_time=False,
				):
		#
		if x0 is not None:
			self.set_x0(x0=x0,
						evaluate_DA=True)
		#
		return self.get_moment_DA(
			D=self.Dk_evaluated,
			Ak=self.Ak_evaluated,
			n=n,
			order=order,
			dimensionless_position=dimensionless_position,
			dimensionless_time=dimensionless_time,
				)



	def get_moment_from_derivatives(self,
			D_derivatives, # numpy array, 
			  # D_derivatives[i] in units of self.D0 / self.L**i
			a_derivatives, # numpy array,
			  # a_derivatives[i] in units of self.L**(1-i)/self.T
			n=0, # n-th moment
			order=4, # order in epsilon for dimensionless (!) moment
			# note that if dimensionless_position == True, then "order"
			# is the order in physical time
			dimensionless_position=False,
			dimensionless_time=False,
				):
		#
		Dk, Ak = self.get_Dk_Ak_from_derivatives(
					D_derivatives=D_derivatives,
					a_derivatives=a_derivatives
						)
		#
		return self.get_moment_DA(Dk=Dk,Ak=Ak,n=n,order=order,
						D_x0=D_derivatives[0],
						dimensionless_position=dimensionless_position,
						dimensionless_time=dimensionless_time)


	def get_moment_DA(self,
				Dk, # numpy arrays, dimensionless quantities
				Ak, # numpy arrays, dimensionless quantities
				n=0, # n-th moment
				order=4, # order in epsilon for dimensionless (!) moment
				# note that if dimensionless_position == True, then "order"
				# is the order in physical time
				dimensionless_position=False,
				dimensionless_time=False,
				D_x0=None,
				):
		#
		# recall that tDL is time, expressed in units of tau_D
		epsilon = sympy.sqrt(2 * self.tDL)
		#
		# construct perturbation series
		moment = 0
		for k in range(order+1):
			moment += epsilon**k * self.moments[n][k]
		#
		# substitute numerical values for Ak, Dk
		moment = self.substitute_Dk_Ak(expr=moment,Dk=Dk,Ak=Ak)
		#
		if not dimensionless_position:
			# convert to length in units of L
			moment = moment * epsilon**n
		#
		if not dimensionless_time:
			# convert to time in units of T
			if D_x0 is None:
				D_x0 = self.D_x0_float
			moment = self.substitute_physical_time(expr=moment,D_x0=D_x0)
			#
		#
		moment = sympy.simplify(moment)
		#
		return sympy.lambdify(self.tDL,moment)


	#############################################################
	# General entropy production power series from coefficients #
	#############################################################

	def get_S_DA(self,
				S_coefficients,
				Dk, # numpy arrays, dimensionless quantities
				Ak, # numpy arrays, dimensionless quantities
				order=4, # order in perturbation series:
				# - if dimensionless_time == True, then "order"
				#   is the order in physical time
				# - if dimensionless_time == False, then "order" 
				#   is the order in epsilon
				dimensionless_time=False,
				D_x0=None,
				include_log_term=False,
				lowest_power=0,
				):
		#
		# recall that tDL is time, expressed in units of tau_D
		epsilon = sympy.sqrt(2 * self.tDL)
		#
		if dimensionless_time:
			range_bound = order+1
		else:
			range_bound = 2*order+1
		#
		# construct perturbation series
		S = 0
		if include_log_term:
			S += sympy.log(epsilon**2)/2
		#
		for k in range(range_bound):
			S += epsilon**(k+lowest_power) * S_coefficients[k]
		#
		# substitute numerical values for Ak, Dk
		S = self.substitute_Dk_Ak(expr=S,Dk=Dk,Ak=Ak)
		#
		if not dimensionless_time:
			# convert to time in units of T
			if D_x0 is None:
				D_x0 = self.D_x0_float
			#
			# the 1/D_x0 is because in dimensionless units, 
			# Sm_dot is units of \tau_D = L**2 / D(x_0).
			# To convert to physical units, we therefore need to
			# multiply by
			# \tau_D / T = D_0/D_x0
			S = self.substitute_physical_time(expr=S * D_x0,
													D_x0=D_x0)
			#
		#
		S = sympy.simplify(S)
		#
		return sympy.lambdify(self.tDL,S)


	##################################
	# Medium entropy production rate #
	##################################

	def get_dS_m_dt(self,
				x0=None,
				order=4,
				dimensionless_time=False,
				):
		#
		if x0 is not None:
			self.set_x0(x0=x0,
						evaluate_DA=True)
		#
		return self.get_S_DA(
			S_coefficients=self.dS_m_dt,
			D=self.Dk_evaluated,
			Ak=self.Ak_evaluated,
			order=order,
			dimensionless_time=dimensionless_time,
				)


	def get_dS_m_dt_from_derivatives(self,
			D_derivatives, # numpy array, 
			  # D_derivatives[i] in units of self.D0 / self.L**i
			a_derivatives, # numpy array,
			  # a_derivatives[i] in units of self.L**(1-i)/self.T
			order=4, # order in epsilon for dimensionless (!) moment
			# note that if dimensionless_position == True, then "order"
			# is the order in physical time
			dimensionless_time=False,
				):
		#
		Dk, Ak = self.get_Dk_Ak_from_derivatives(
					D_derivatives=D_derivatives,
					a_derivatives=a_derivatives
						)
		#
		return self.get_S_DA(
						S_coefficients=self.dS_m_dt,
						Dk=Dk,Ak=Ak,order=order,
						D_x0=D_derivatives[0],
						dimensionless_time=dimensionless_time)


	#################################
	# Total entropy production rate #
	#################################

	def get_dS_tot_dt(self,
				x0=None,
				order=4,
				dimensionless_time=False,
				):
		#
		if x0 is not None:
			self.set_x0(x0=x0,
						evaluate_DA=True)
		#
		return self.get_S_DA(
			S_coefficients=self.dS_tot_dt,
			D=self.Dk_evaluated,
			Ak=self.Ak_evaluated,
			order=order,
			dimensionless_time=dimensionless_time,
			lowest_power=-2,
				)


	def get_dS_tot_dt_from_derivatives(self,
			D_derivatives, # numpy array, 
			  # D_derivatives[i] in units of self.D0 / self.L**i
			a_derivatives, # numpy array,
			  # a_derivatives[i] in units of self.L**(1-i)/self.T
			order=4, # order in epsilon for dimensionless (!) moment
			# note that if dimensionless_position == True, then "order"
			# is the order in physical time
			dimensionless_time=False,
				):
		#
		Dk, Ak = self.get_Dk_Ak_from_derivatives(
					D_derivatives=D_derivatives,
					a_derivatives=a_derivatives
						)
		#
		return self.get_S_DA(
						S_coefficients=self.dS_tot_dt,
						Dk=Dk,Ak=Ak,order=order,
						D_x0=D_derivatives[0],
						dimensionless_time=dimensionless_time,
						lowest_power=-2,
						)


	#################
	# Gibbs entropy #
	#################

	def get_S_Gibbs(self,
				x0=None,
				order=4,
				dimensionless_time=False,
				):
		#
		if x0 is not None:
			self.set_x0(x0=x0,
						evaluate_DA=True)
		#
		return self.get_S_DA(
			S_coefficients=self.S_Gibbs,
			D=self.Dk_evaluated,
			Ak=self.Ak_evaluated,
			order=order,
			dimensionless_time=dimensionless_time,
			include_log_term=True,
				)


	def get_S_Gibbs_from_derivatives(self,
			D_derivatives, # numpy array, 
			  # D_derivatives[i] in units of self.D0 / self.L**i
			a_derivatives, # numpy array,
			  # a_derivatives[i] in units of self.L**(1-i)/self.T
			order=4, # order in epsilon for dimensionless (!) moment
			# note that if dimensionless_position == True, then "order"
			# is the order in physical time
			dimensionless_time=False,
				):
		#
		Dk, Ak = self.get_Dk_Ak_from_derivatives(
					D_derivatives=D_derivatives,
					a_derivatives=a_derivatives,
						)
		#
		return self.get_S_DA(
						S_coefficients=self.S_Gibbs,
						Dk=Dk,Ak=Ak,order=order,
						D_x0=D_derivatives[0],
						dimensionless_time=dimensionless_time,
						include_log_term=True,
						)


	##################################
	# Some internal helper functions #
	##################################

	def __get_N_max(self):
		''' Return index of largest appearing Dk, Ak in the pre-calculated 
		symbolic definitions

		We use this function to determine how many symbolic variables for 
		the D_k, A_k we need to define upon instantiation of the class

		''';
		#
		# this should always be a complete list of all the precalculated
		# symbolic expressions that are used in this module, and imported
		# at line 7 and the following lines
		list_of_symbolic_definitions = [
						npp,
						ppp,
						moments,
						dS_m_dt,
						dS_tot_dt,
						S_Gibbs,
						]
		#
		N_max = 0
		for namespace in list_of_symbolic_definitions:
			N_max = max(N_max, namespace.N_max)
		return N_max

	def __eval_D_x0(self,x0):
		#
		x0 = self.check_dimensions(quantity=x0,
									dimensions=self.L,
									return_with_dimensions=False)
		#
		self.D_x0_float = self.D_lambda(x0)
		self.D_x0 = self.D_x0_float * self.D0

	def __eval_Dk(self,x0):
		''' Evaluates lambda functions for D_n for given position x_0

		This method defines 1D numpy.arrays 
			self.Dk_evaluated
		of length self.N_Max, such that 
			self.Dk_evalauted[n] = D_n
		with D_n defined in Eq. (29) of Ref. [1]

		Args:
			x0 (numpy.float): Position in units of L
		Returns: None
		'''
		#
		if not self.D_functions_constructed:
			# if self.Dk_lambda has not been defined yet, try to define it
			self.construct_D_functions()
		#
		self.Dk_evaluated = np.zeros(self.N_max,dtype=float)
		#
		for i, function in enumerate(self.Dk_lambda):
			self.Dk_evaluated[i] = function(x0)


	def __eval_Ak(self,x0):
		''' Evaluates lambda functions for A_n for given position x_0

		This method defines 1D numpy.arrays 
			self.Ak_evaluated
		of length self.N_Max, such that 
			self.Ak_evalauted[n] = A_n
		with A_n defined in Eq. (28) of Ref. [1]

		Args:
			x0 (numpy.float): Position in units of L
		Returns: None
		'''
		#
		if not self.A_functions_constructed:
			# if self.Ak_lambda has not been defined yet, try to define it
			self.construct_A_functions()
		#
		self.Ak_evaluated = np.zeros(self.N_max,dtype=float)
		#
		for i,current_function in enumerate(self.Ak_lambda):
			self.Ak_evaluated[i] = current_function(x0)



	def __substitute_local_variables(self,expression,remote):
		'''Substitute local symbolic variables into imported expressions

		All perturbative analytical results are stored in separate python
		files (see subdirectory sympy_definitions/), and are imported at 
		the beginning of this file (starting in line 7).
		Upon importing, each analytical expression comes with its own symbolic
		variables for space, time, etc

		This function takes an expression and substitutes the remote symbolic
		variables with the local symbolic variables used in this class.

		This function is class private as normally there is no reason for the
		user to call it.

		Args:
			expression (sympy expression): sympy expression that uses symbolic
										   variables from the remote namespace
			remote (sympy expression): remote namespace

		Returns:
			sympy expression: the input expression, but with symbolic variables
							  from the self namespace (i.e. from this class)
		'''
		#
		if not isinstance(expression,sympy.Expr):
			return expression
		#
		try:
			expression = expression.subs(remote.xDL,self.xDL)
		except AttributeError:
			# this means that in the file of the expression, the symbol
			# xDL has not been defined. This is not a problem, so we just
			# continue
			pass
		#
		try:
			expression = expression.subs(remote.epsilon,self.epsilon)
		except AttributeError:
			# this means that in the file of the expression, the symbol
			# epsilon has not been defined. This is not a problem, so we just
			# continue
			pass
		#
		for j in range(self.N_max):
			#
			try:
				expression = expression.subs(remote.Dk[j],self.Dk[j])
			except IndexError:
				# this means that the expression doesn't contain the Dk
				# up self.N_max. This by itself is not a problem, so we
				# just continue
				pass
			except AttributeError:
				# this means that in the file of the expression, the list
				# Dk has not been defined. This is not a problem, so we just
				# continue
				pass
			#
			try:
				expression = expression.subs(remote.Ak[j],self.Ak[j])
			except IndexError:
				# this means that the expression doesn't contain the Ak
				# up self.N_max. This by itself is not a problem, so we
				# just continue
				pass
			except AttributeError:
				# this means that in the file of the expression, the list
				# Ak has not been defined. This is not a problem, so we just
				# continue
				pass
		#
		return expression
