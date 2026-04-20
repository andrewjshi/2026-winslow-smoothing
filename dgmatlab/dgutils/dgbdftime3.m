function u=dgbdftime(u,msh,data,phys,jacname,tpars,spars,solver,pltpars,filename)
%DGBDFTIME  BDF Time Stepping, upto 3rd order
%
%    Solve nonlinear problem with damped Newton iterations.
%
%    Syntax: u=dgbdftime(u,msh,data,phys,jacname,tpars,spars,solver,pltpars,filename)
%
%    u:                Initial guesses (input)
%    u:                Solutions (output)
%    msh, data, phys:  DG problem
%    jacname:          Jacobian/Residual function
%    tpars:            = [dt, nsteps, pltfreq, maxorder]
%                dt:     Timestep
%            nsteps:     Number of timesteps
%           pltfreq:     Plotting frequency (inf for no plotting)
%          maxorder:     Maximum order of timestepping
%    spars:            = [newtontol, newtonmaxiter, linsystol, reportlevel]
%         newtontol:     Residual tolerance in Newton iterations (1e-8)
%     newtonmaxiter:     Maximum number of Newton steps           (20)
%         linsystol:     Residual tolerance in iterative solver  (1e-5)
%     linsysmaxiter:     Maximum number of Krylov iterations     (1000)
%       reportlevel:     Level of diagnostic output in linear solver (0)
%    solver:           Solver command
%                      = 'direct'     : Direct solver (default)
%                        'gmres j'    : GMRES, block Jacobi precond.
%                        'gmres i'    : GMRES, block ILU(0) precond.
%                        'gmres 0dpi' : GMRES, ILU(0)-p0 precond.
%                        'gmres'      :       -    "    -
%
%    Optional plotting parameters:
%
%    pltpars:          Parameters to dgplot (cell array)
%                      Empty for no plotting
%    filename:         Filename to save solution (for recovery)

if nargin<7 | isempty(spars), spars=[]; end
if nargin<8 | isempty(solver), solver='direct'; end
if nargin<9, pltpars={}; end
if nargin<10, filename=''; end

if length(spars)<1, spars=[spars,1e-8]; end
if length(spars)<2, spars=[spars,20]; end
if length(spars)<3, spars=[spars,1e-5]; end
if length(spars)<4, spars=[spars,1000]; end
if length(spars)<5, spars=[spars,0]; end

if length(tpars)<3,tpars=[tpars,inf]; end
if length(tpars)<4, tpars=[tpars, 1]; end



ix=find(solver==' ');
if isempty(ix)
    solvecommand=solver;
    precond='';
else
    solvecommand=solver(1:ix-1);
    precond=solver(ix+1:end);
end

dt=tpars(1);
nsteps=tpars(2);
maxorder = tpars(4);
pltfreq=tpars(3);

iorder = size(u,4); % initial order
order = iorder;     % current order

if ~isfield(phys,'time'), phys.time=0.0; end

for ii=1:nsteps
    newtontol=spars(1);
    maxiter=spars(2);

    dim=size(msh.s,2)-1;

    phys.time=phys.time+dt;

    R=feval(jacname,u,msh,data,phys);
    normR=dt*max(abs(R(:)));


    if ii==1, N=size(R,2); M=mkmassmatrix(msh,data,N,true);
    end

    % set uhis initially from input solns
    if ( order == iorder && ii < 2 ), 
        uhis = u; 
        clear u, 
        u = uhis(:,:,:,end);  
    end

   
    
    switch order
        case 1

            fprintf('Timestep %4d,\t BDF order %4d\n',ii,order);

            clear DA OA
            R=feval(jacname,u,msh,data,phys);
            normR=dt*max(abs(R(:)));
            
            u0=u;
            uvar=u0; coef=1;

            fprintf('  Newton iter %2d, Res = %18.14f',0,normR);

            [u,normR]=newtonsolve(jacname,u,uvar,coef,dt,msh,data,phys,...
                                       precond,spars,solvecommand,solver,M);
           
            if order < maxorder, 
                order = order + 1; 
                uhis(:,:,:,order) = u; 
                uhis(:,:,:,1) = u0;
            else
                uhis(:,:,:,order) = u;
            end

            
        case 2

            fprintf('Timestep %4d,\t BDF order %4d\n',ii,order);
            u1 = uhis(:,:,:,1);
            u0 = uhis(:,:,:,2);

            clear DA OA
            R=feval(jacname,u,msh,data,phys);
            normR=dt*max(abs(R(:)));

            uvar=(4*u0-u1)/3; coef=2/3;
            
            
            fprintf('  Newton iter %2d, Res = %18.14f',0,normR);

            [u,normR]=newtonsolve(jacname,u,uvar,coef,dt,msh,data,phys,...
                                        precond,spars,solvecommand,solver,M);
           
            % update solutions
            if order < maxorder,
                order = order + 1;
                uhis(:,:,:,order) = u; 
                uhis(:,:,:,2) = u0;
                uhis(:,:,:,1) = u1;
            else 
                uhis(:,:,:,2) = u;
                uhis(:,:,:,1) = u0;
            end

           
            
        case 3
            fprintf('Timestep %4d,\t BDF order %4d\n',ii,order);
            u2 = uhis(:,:,:,1);
            u1 = uhis(:,:,:,2);
            u0 = uhis(:,:,:,3);
      
            uvar=(18*u0-9*u1+2*u2)/11; coef=6/11;

            fprintf('  Newton iter %2d, Res = %18.14f',0,normR);

            [u,normR]=newtonsolve(jacname,u,uvar,coef,dt,msh,data,phys,...
                                        precond,spars,solvecommand,solver,M);
            
            % update solutions
            uhis(:,:,:,3) = u;
            uhis(:,:,:,2) = u0;
            uhis(:,:,:,1) = u1;

            
            
    
    end
    fprintf('\n');
    if normR>=newtontol
        warning('No convergence.');
    end

    if ~isempty(pltpars) & mod(ii,pltfreq)==0
        dgplot(msh,u,pltpars{:});
    end
    if ~isempty(filename)
        save(filename,'u');
    end
    
   
end
u = uhis;


function [u,normR]=newtonsolve(jacname,u,uvar,coef,dt,msh,data,phys,...
                                    precond,spars,solvecommand,solver,M)

newtontol=spars(1);
maxiter=spars(2);

for jj=1:maxiter
    if jj>1, clear DA OA, end
    [R,DA,OA]=feval(jacname,u,msh,data,phys);
    R=dgmassinv(u-uvar,msh,data,1)-coef*dt*R;
    DA=M-coef*dt*DA;
    OA=-coef*dt*OA;

    [du,outpars]=dglinsys(msh,data,DA,OA,R,precond,spars(3:5),solvecommand);
    if strcmp(solver,'direct')
        fprintf('\n');
    else
        fprintf(', Krylov iters = %4d, res = %5.2g\n',outpars(2),outpars(1));
    end
    if any(isnan(du(:))), error('NaN in solution.'); end
    u=u-du;


    clear DA OA
    R=feval(jacname,u,msh,data,phys);
    R=dgmassinv(u-uvar,msh,data,1)-coef*dt*R;
    normR=max(abs(R(:)));
    fprintf('  Newton iter %2d, Res = %18.14f',jj,normR);
    if normR<newtontol, break; end
end






