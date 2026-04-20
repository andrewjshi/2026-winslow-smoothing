function u = dglinsys(msh,data,DA,OA,b,precond,pars,solver)
%DGLINSYS Solve linear system of equations, DG block format
%
%   U = DGLINSYS(MSH,DATA,DA,OA,B,PRECOND,PARS,SOLVER) solves
%   the linear system AU=B using SOLVER preconditioned with
%   PRECOND. The matrix A is represented in block format by DA,OA.
%   PARS = [TOL,REPORTLEVEL], where TOL is residual tolerance and
%   the integer REPORTLEVEL >= 0 is the amount of diagnostic output.
%
%   Defaults: MSH, DATA, DA, OA, B must be given
%             PRECOND = '0dpi'
%             PARS    = [TOL, MAXITER, REPORTLEVEL] = [1e-6, 1000, 2]
%             SOLVER  = 'gmres'
%
%   Common preconditioners:
%             'j'     : Block Jacobi
%             'i'     : Block ILU(0)
%             'j0dpj' : Jacobi/p=0 correction
%             '0dpi'  : ILU(0)/p=0 correction
