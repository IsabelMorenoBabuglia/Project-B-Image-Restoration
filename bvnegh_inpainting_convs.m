function    bvnegh_inpainting_convs

% TV-H{-1} inpainting based on:
% [1] M. Burger, L. He, C.-B. Schönlieb, Cahn-Hilliard inpainting and a
% generalization for grayvalue images, SIAM J. Imaging Sci. Volume 2, 
% Issue 4, pp. 1129-1167 (2009).
% 
% It was discretized based on a convexity splitting
% method
% [2] C.-B. Schönlieb, A. Bertozzi, Unconditionally stable schemes for
% higher order inpainting, Communications in Mathematical Sciences,
% Volume 9, Issue 2, pp. 413-457 (2011).

% namely:
% E_1  = \int_{\Omega} |\nabla u| dx ,
% E_11 = \int_{\Omega}  C_1/2 |\nabla u|^2 dx, E_12 = \int_{\Omega} - |\nabla u| + C_1/2 |\nabla u|^2 dx

% E_2 = \lambda \int_{\Omega\D} (f-u)^2 dx,
% E_21 = \int_{\Omega\D} C_2/2 |u|^2 dx, E_22 = \int_{\Omega\D} -\lambda (f-u)^2 + C_2/2|u|^2 dx

%December 2011
%Carola Schönlieb

%Reading in and preparing image data
load strip1ades.mat
g=strip1ades;
[m,n]=size(g);
g=double(g);
g = g./max(max(g));
clims=[0 1];
figure(1)
imagesc(g,clims); axis image; axis off; colormap(gray);

%%%%%%%%%%%%%%%%%%%% Definition of PARAMETERS:%%%%%%%%%%%%%%%%%%%%%%%
h1=1;
h2=1;
dt=100;
epsilon= 0.01;
ep2 = epsilon^2;
lpower=1;
lambda0 = 10^lpower;
lambda = lambda.*lambda0;
c1=1/epsilon;
c2=lambda0;
Itermax = 2000;

% Diagonalize the Laplace Operator by: Lu + uL => D QuQ + QuQ D, where Q is nonsingular, the matrix of eigenvectors of L
% and D is a diagonal matrix.
% We have to compute QuQ. This we can do in a fast way by using the fft-transform:

% Initialization of u and its Fourier transform:
u=g;
ubar=fft2(u);
lu0bar=fft2(lambda.*u);
lubar=lu0bar;
curv=zeros(m,n);

%Definition of the Laplace operator with periodic boundary conditions:
L1=1/(h1^2)*(diag(-2*ones(m,1)) + diag(ones(m-1,1),1) + diag(ones(m-1,1),-1)...
    + diag(ones(1,1),m-1)+ diag(ones(1,1),-m+1));
L2=1/(h2^2)*(diag(-2*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1)...
    + diag(ones(1,1),n-1)+ diag(ones(1,1),-n+1));


% Computation of the eigenvalues of L1:
for j=1:m
    eigv1(j) = 2*(cos(2*(j-1)*pi/m)-1);
    Lambda1(j,j)=eigv1(j);
end 

% Computation of the eigenvalues of L2:
for j=1:n
    eigv2(j) = 2*(cos(2*(j-1)*pi/n)-1);
    Lambda2(j,j)=eigv2(j);
end    

Denominator=1/h1^2.*Lambda1*ones(m,n) + 1/h2^2.*ones(m,n)*Lambda2;

% Now we can write the above equation in  much simpler way and compute the solution ubar

it=0;
err=1;
while it < Itermax   
    
    % Computation of the tv-seminorm:
    % estimate derivatives
    ux = (u(:,[2:n n])-u(:,[1 1:n-1]))/2;
	uy = (u([2:m m],:)-u([1 1:m-1],:))/2;
	uxx = u(:,[2:n n])+u(:,[1 1:n-1])-2*u;
	uyy = u([2:m m],:)+u([1 1:m-1],:)-2*u;
	Dp = u([2:m m],[2:n n])+u([1 1:m-1],[1 1:n-1]);
	Dm = u([1 1:m-1],[2:n n])+u([2:m m],[1 1:n-1]);
	uxy = (Dp-Dm)/4;
    % compute flow
    Num = uxx.*(ep2+uy.^2)-2*ux.*uy.*uxy+uyy.*(ep2+ux.^2);
    Den = (ep2+ux.^2+uy.^2).^(3/2);
    curv = Num./Den;
    lapcurvbar = Denominator.*fft2(curv);
    
    ubar = ((1+c2*dt+dt*c1*Denominator.^2).*ubar - dt*lapcurvbar ...
        + dt*(lu0bar-lubar))./(1+dt*c2+dt*c1*Denominator.^2);
    u = real(ifft2(ubar));
    it=it+1;
    figure(2);
    imagesc(u,clims); axis image; axis off; colormap(gray);
    title(['Inpainted image as solution to the fourth0-order total '...
        'variation flow after ' num2str(it) ' timesteps']);
    pause(0.01)
    lubar=fft2(lambda.*u);
end