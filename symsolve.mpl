# Last updated: July 31 2017
# Author: Caroline Lin
# Written for Maple 2016, on Ubuntu Linux 16.04 LTS.
#
# Typical usage
# 
# Comments
# * for approximate case, inequalities must be specified in terms of the approximated functions
#   e.g. f(u) <> 3 should be replaced by f(v[0][0]) <> 3
#   otherwise rifsimp will error later on. 


with(PDEtools):
with(DEtools):

printCasenum := proc(rifcases)
	return rifcases[casecount];
end proc;

# print the output of rifsimp with cases in an easy-to-read manner.
printCases := proc(rifcases, limit := NULL, verbose := true, lowlim := 1)
	local i, uplim;
	uplim := min(limit, rifcases[casecount]);
	for i from lowlim to uplim do 
		print("=================================== CASE", i, "=========================================");
		
		if verbose = true then
			printEqns(rifcases[i][Solved]);
		end if;

		print(rifcases[i][Pivots]); 
		print(DEtools[initialdata](rifcases[i][Solved]));
	end do;
end proc;

# print the output of a list of equations.
printEqns := proc(equations)
	local i;
	for i from 1 to nops(equations) do
		print(equations[i]);
	end do;
end proc;


# top-level single PDE. system is determining system, rifresults are the result of rifsimp.
symSolve := proc(x, y, eqn, ineqs := [], p := 1, epsilon := NULL, cases := false)
	local system, tosolve, rifresults, i;
	if epsilon <> false then
		printf("CHECK: Approximate case \n");
		system, tosolve := op(approximate(x, y, eqn, p, epsilon));
	else
		printf("CHECK: Exact case \n");
		system, tosolve := eqn, y;	
	end if;
	
	if epsilon <> false then print("APPROXIMATE: System has been approximated to order ", p); end if;
	print("SOLVING FOR:", tosolve);
	system := determine(x, tosolve, system); 
	rifresults := _eliminate([op(x), op(tosolve)], system[2], [op(system[1]), op(ineqs)], cases); 
	return system, rifresults;
end proc;

# top-level function for systems of PDEs.
systemOverhead := proc()
	printf("Function not yet implemented");
end proc;


# replace PDE with power-series approximation in epsilon.
# The unknown functions in list 'y' are replaced with coefficient functions.
approximate := proc(x::list, y::list, eqn, p, epsilon)
	local modified, i, pol;
	global v; # v is global because I'm paranoid about local values of v interfering with
	unassign(v); # the return of this function
	modified := eqn;
	
	if not has(eqn, y) and not has(eqn, x) then
		WARNING("Supplied equation %1 does not contain the specified variables", eqn);
	end if;
	
	# Substitute order-p polynomial approximations into 'modified'.
	for i from 1 to nops(y) do
		pol[i] := add(v[i - 1][j](op(x)) * epsilon^j, j=0..p); # v is just an arbitrary table of functions.
		modified := subs(y[i](op(x)) = pol[i], modified);
	end do;

	modified := lhs(modified) - rhs(modified); # conv. eqn to expr
	modified := taylor(modified, epsilon, p + 1); # Make taylor series in epsilon. Groups terms by power of epsilon.
	modified := convert(modified, polynom); # delete the Order term by converting type Series -> Polynomial
	
	return [extractSystem(modified, epsilon), componentFunctions(v, nops(y) - 1, p)];
end proc;

# (hidden) helper method
componentFunctions := proc(v, row, col)
	return [seq(seq(v[i][j], i = 0 .. row), j = 0 .. col)];
end proc;

# (hidden) helper method
# The coefficients of a perturbed PDE are also PDEs. Extract the set (system) of PDEs given by the coefficients.
extractSystem := proc(approx_pert, epsilon)
	local system;
	system := [coeffs(approx_pert, epsilon)];
	system := [seq(system[i] = 0, i = 1 .. nops(system))]; # turn exprs into eqns.
	return system;
end proc;


# Call DeterminingSystem (linear conditions for symmetries of PDE system).
determine := proc(x, y, system)
	local xis, etas, generators;
	global xi, eta; 
	unassign(xi); unassign(eta); # reserve xi and eta.
	
	xis := [seq(xi[i](op(x), op(y)), i = 1 .. nops(x))];
	etas := [seq(eta[i](op(x), op(y)), i = 1 .. nops(y))];
	generators := [op(xis), op(etas)];	
	userinfo(2, `GENERATORS:`, generators);
	return [DeterminingPDE(system, generators, integrabilityconditions = false), generators];
end proc;

# Call rifsimp. Wrapper deals with multiple cases.
_eliminate := proc(x, y, system, cases := false)
	local rifresults;
	# syntax: rifsimp( system, dependent varables, options.)
	if cases = true then 
		rifresults := DEtools[rifsimp](system, y, casesplit);
	else
		rifresults := table([1 = DEtools[rifsimp](system, y), casecount = 1]);
	end if;
	return rifresults;
end proc;

# [Not in use] used to be called by _eliminate, but pdsolve errors out too much
simplesolve := proc(x, y, system)
	local solutions, degfreedom, _initialdata;
	solutions := NULL; degfreedom := NULL;
	 	
	_initialdata := DEtools[initialdata](system);	
	try	
		solutions := pdsolve(system, y); # syntax: pdsolve(system, y)	
	catch:
		print("************* error **************");
		print("_ELIMINATE: error occurred while trying to pdsolve");
		print(lastexception);
		print("************* end error *************");
	end try;
	return [solutions, _initialdata];
end proc:


