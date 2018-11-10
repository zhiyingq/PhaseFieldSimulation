%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase Field Finite Difference %
%           Code For            %
%   Dendritic Crystallization   %
%   Optimized For Matlab        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%== get initial wall time: (format long:显示15位双精度)
time0=clock();
format long;

%-- Simulation cell parameters:
Nx=300;     %x格点数
Ny=300;     %y格点数
NxNy=Nx*Ny;

dx=0.03;    %grid spacing between two grid points in the x-direction
dy=0.03;    %grid spacing between two grid points in the x-direction  

%-- Time integration parameters:

nstep=4000;     % time integration seps
nprint=50;      % Output frequency to write the results to file
dtime=1.e-4;    % Time increment for numerical integration

%-- Material specific parameters:

tau=0.0003;      % τ
epsilonb=0.01;   % ε拔
mu=1.0;          % μ
kappa=1.8;       % κ
delta=0.02;      % δ
aniso=6.0;       % j
alpha=0.9;       % α
gamma=10.0;      % γ
teq=1.0;         % Teq
theta0=0.2;      % θ0
seed=5.0;        % The size of the initial seed (In grid numbers) 
                 % See function nucleus
pix=4.0*atan(1.0);  % The value of pi

%--- Initialize and introduce initial nuclei;

[phi,tempr]=nucleus(Nx,Ny,seed);
% Initialize the phi and temper arrays and introduce the seed in center
%--- Laplacian templet

[laplacian]=laplacian(Nx,Ny,dx,dy);
%Calculate the finite difference template for the laplacians

%-- Evolution, calculate Eqs4.51 and Eqs4.52
%--

for istep=1:nstep
    
    phiold=phi; %首先定义了phi_old和phi是相等的
    
    %--
    % calculate the laplacians and epsilon:
    %--
    
    phi2=reshape(phi',NxNy,1); %将phi(Nx,Ny)转换成一维数组phi2(Nx*Ny)，phi'为phi的转置，只
                               %只有转置后才能正确地按照第一行、第二行的顺序进行reshape
    lap_phi2=laplacian*phi2;   %计算phi2的拉普拉斯值，?2φ
    
    [lap_phi]=vec2matx(lap_phi2,Nx);    %将lap_phi2(Nx*Ny）转化到二维矩阵lap_phi(Nx,Ny)
    
    %--
    
    tempx=reshape(tempr',NxNy,1); %将tempr(Nx,Ny)转换为一维数组tempx(NxNy)
    
    lap_tempx=laplacian*tempx; %计算?2T
    
    [lap_tempr]=vec2matx(lap_tempx,Nx); %将lap_tempx(Nx*Ny)转换为二维矩阵lap_tempr(Nx,Ny)
    
    %--gradients of phi:
    
    [phidy,phidx]=gradient_mat(phi,Nx,Ny,dx,dy); %计算Cartesian derivatives of phi,?φ
    %一维方向上的梯度，但加上了边界条件以及考虑了dx和dy
    %calculate angle:
    
    theta=atan2(phidy,phidx); %calculate theta:θ=tan-1(φy/φx)
    %这里的atan2是定义在(-π,π)上的，这里的theta是一个数组
    %--epsilon and its derivative:
    
    epsilon=epsilonb*(1.0+delta*cos(aniso*(theta-theta0))); %calculate ε,它是一个数组
    
    epsilon_deriv=-epsilonb*aniso*delta*sin(aniso.*(theta-theta0));
    %calculate ?ε/?θ  .*指的是两个矩阵的对应元素相乘，aij*bij,而不是矩阵乘法，它是一个数组
    %--- first term of Eq.4.51 计算第一项
    
    dummyx=epsilon.*epsilon_deriv.*phidx;
    
    [term1,dummy]=gradient_mat(dummyx,Nx,Ny,dx,dy);
    
    %--- second term of Eq.4.51 计算第二项
    
    dummyy=-epsilon.*epsilon_deriv.*phidy;
    
    [dummy,term2]=gradient_mat(dummyy,Nx,Ny,dx,dy);
    
    %--- factor m:
    
    m=(alpha/pix)*atan(gamma*(teq-tempr));
    
    %-- Time integration:
    
    phi=phi+(dtime/tau)*(term1+term2+epsilon.^2.*lap_phi+...
        phiold.*(1.0-phiold).*(phiold-0.5+m));
    
    %-- evolve temperature:
    
    tempr=tempr+dtime*lap_tempr+kappa*(phi-phiold);
    
    
    %---- print results
    
    
    if(mod(istep,nprint)==0)
        
        fprintf('done step: %5d\n',istep);
        
        fname1=sprintf('time_%d_phi.txt',istep);
        out1=fopen(fname1,'w');
        
        for i=1:Nx
            for j=1:Ny
             
                fprintf(out1,'%14.6e\n',phi(i,j));
            end
        end
        
        fname2=sprintf('time_%d_tempr.txt',istep);
        out2=fopen(fname2,'w');
        
        for i=1:Nx
            for j=1:Ny
             
                fprintf(out2,'%14.6e\n',tempr(i,j));
            end
        end
      
        
        %write_vtk_grid_values(Nx,Ny,dx,dy,istep,phi,tempr);
        
    end
end

%calculate compute time:

compute_time=etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);

    
    
    
    
    
    
    
    
    


