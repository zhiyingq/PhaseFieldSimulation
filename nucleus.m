%This function introduces initial solid nuclei in the 
%center of the simulation cell. The size of the nuclei
%is in terms of grid numbers
%其实就是一个初始化函数，Φ和T是两个二维数组，初始值都为0
%以网格中心为原点，半径小于seed的的都认为是种子，这时Φ=1
%输出Φ和T两个二维数组

function[phi,tempr]=nucleus(Nx,Ny,seed)

format long;

phi=zeros(Nx,Ny);
tempr=zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        if((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2)<seed)
            phi(i,j)=1.0;
        end
    end
end

end
