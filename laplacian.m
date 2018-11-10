%This function generates the Laplae operator by using five-point
%stencil with periodic boundaries. It is used in optimized Matlab
%得到的是一个nx平方*ny平方的矩阵，例如，本来是10X10的矩阵，变成了100X100
%用这个grad矩阵乘以一个(100X1)的矩阵，即包含了每个点的列，则得到了每一个点的梯度值
%该梯度值等于-4*x+(上+下+左+右），采用循环边界，最上的上=最下，最左的左=右，以此类推

function[grad]=laplacian(nx,ny,dx,dy)

format long;

nxny=nx*ny; %total number of grid points in the simulation cell

r=zeros(1,nx);
r(1:2)=[2,-1];
T=toeplitz(r);

E=speye(nx); %创造一个nx*nx的对角线为1的稀疏矩阵的

grad=-(kron(T,E)+kron(E,T));
%Kron即为Kronecker积，X,Y的Kron积为：
%X11*Y X12*Y .....X1N*Y
%X21*Y X22*Y .....X2N*Y
%......
%XN1*Y XN2*Y .....XNN*Y

%-- for periodic boundaries

for i=1:nx
    ii=(i-1)*nx+1;
    jj=ii+nx-1;
    grad(ii,jj)=1.0;
    grad(jj,ii)=1.0;
    
    kk=nxny-nx+i;
    grad(i,kk)=1.0;
    grad(kk,i)=1.0;
end

grad=grad/(dx*dy);

end
