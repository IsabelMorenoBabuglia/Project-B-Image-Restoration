function heat_equation
%%%% Inpainting using the heat equation %%%%%%
%December 2011
%Carola Schönlieb

%Reading in and preparing image data
load strip1ades.mat
g=strip1ades;
[m,n]=size(g);
g=double(g);
g = g./max(max(g));
clims=[0 1];
%Parameters
h1=1;
h2=1;
T=1000;
dt=0.1;
lambda0=10;
lambda = lambda.*lambda0;
figure(1); imagesc(g, clims); colormap(gray)

% laplacian = L1*u+u*L2
L1=(1/(h1^2))*(diag(-2*ones(m,1)) + diag(ones(m-1,1),1) + diag(ones(m-1,1),-1));
L1(1,1)=-1/h1^2; L1(m,m)=-1/h1^2;
L2=(1/(h2^2))*(diag(-2*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
L2(1,1)=-1/h2^2; L2(n,n)=-1/h2^2;

%intialization of u:
u=g;
% Explicit time stepping scheme for the heat equation
for t=1:T
    u=u+dt*(L1*u+u*L2+lambda.*(g-u));
    figure(2);
    imagesc(u, clims), axis equal; axis off; colormap(gray)
    title(['Inpainted image as solution to heat equation after '...
        num2str(t) ' timesteps']);
    pause(0.01)
end