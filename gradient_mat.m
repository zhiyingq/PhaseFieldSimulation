%This function calculates the one-dimensinal gradients for x and y
%directions of a given matrix.

function[matdx,matdy]=gradient_mat(matx,Nx,Ny,dx,dy)

format long;

% the built-in function "gradient" is used in this function

[matdx,matdy]=gradient(matx);

%-- for periodic boundaries:

matdx(1:Nx,1)=0.5*(matx(1:Nx,2)-matx(1:Nx,Nx));
matdx(1:Nx,Nx)=0.5*(matx(1:Nx,1)-matx(1:Nx,Nx-1));

matdy(1,1:Ny)=0.5*(matx(2,1:Ny)-matx(Ny,1:Ny));
matdy(Ny,1:Ny)=0.5*(matx(1,1:Ny)-matx(Ny-1,1:Ny));

matdx=2.0*matdx/dx;
matdy=2.0*matdy/dy;

end

