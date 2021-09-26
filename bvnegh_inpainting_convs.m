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
load strip1ades1.mat
g=strip1ades1;
[m,n]=size(g);
g=double(g);
g = g./max(max(g));
clims=[0 1];
figure(1)
imagesc(g,clims); axis image; axis off; colormap(gray);

%%%%%%%%%%%%%%%%%%%% Definition of PARAMETERS:%%%%%%%%%%%%%%%%%%%%%%%
h1=1;                     %h1 is dx
h2=1;                     %h2 is dy
dt=100;
epsilon= 0.01;
ep2 = epsilon^2;
lpower=1;
lambda0 = 10^lpower;
lambda = lambda.*lambda0;   %lambda is defined in the file striplades1.mat , in inpainting area lambda=0, other area lambda=constant
c1=1/epsilon;
c2=lambda0;
Itermax = 2000;

% Diagonalize the Laplace Operator by: Lu + uL => D QuQ + QuQ D, where Q is nonsingular, the matrix of eigenvectors of L
% and D is a diagonal matrix.
% We have to compute QuQ. This we can do in a fast way by using the fft-transform:

%here the operator is split into two parts because u is a rectangle matrix
%and you can't use a square matrix on a rectangle matrix.
%to make sure the operator makes effect on every elements of the matrix, it
%has to be split into two parts

% Initialization of u and its Fourier transform:
%here we use this transform because it can change PDE into algebraic equation
u=g;
ubar=fft2(u);
lu0bar=fft2(lambda.*u);
lubar=lu0bar;
curv=zeros(m,n);

%Definition of the Laplace operator with periodic boundary conditions:
%matrix method, has been written on the note
L1=1/(h1^2)*(diag(-2*ones(m,1)) + diag(ones(m-1,1),1) + diag(ones(m-1,1),-1)...
    + diag(ones(1,1),m-1)+ diag(ones(1,1),-m+1));
L2=1/(h2^2)*(diag(-2*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1)...
    + diag(ones(1,1),n-1)+ diag(ones(1,1),-n+1));


% Computation of the eigenvalues of L1:
% the method to calculate the eigenvalue of this matrix comes from 
% https://click.endnote.com/viewer?doi=10.1007%2F3-540-44935-3_16&token=WzM0MTEyNTIsIjEwLjEwMDcvMy01NDAtNDQ5MzUtM18xNiJd.wtVs4CZWek24z-3PjCJy15jyvFQ
for j=1:m
    eigv1(j) = 2*(cos(2*(j-1)*pi/m)-1);
    Lambda1(j,j)=eigv1(j);
end 

% Computation of the eigenvalues of L2:
for j=1:n
    eigv2(j) = 2*(cos(2*(j-1)*pi/n)-1);
    Lambda2(j,j)=eigv2(j);
end    

%here Denominator is a diagonal matrix, you can regard it as laplace
%operator in the fourier space
Denominator=1/h1^2.*Lambda1*ones(m,n) + 1/h2^2.*ones(m,n)*Lambda2;

% Now we can write the above equation in  much simpler way and compute the solution ubar

it=0;
err=1;
while it < Itermax   
    
    % Computation of the tv-seminorm:
    % estimate derivatives
    ux = (u(:,[2:n n])-u(:,[1 1:n-1]))/2;    %for example, if u=[1 2 3], u(:,[2:3 3])=[2 3 3], then you can see why it can express ux
	uy = (u([2:m m],:)-u([1 1:m-1],:))/2;
	uxx = u(:,[2:n n])+u(:,[1 1:n-1])-2*u;
	uyy = u([2:m m],:)+u([1 1:m-1],:)-2*u;
	Dp = u([2:m m],[2:n n])+u([1 1:m-1],[1 1:n-1]);
	Dm = u([1 1:m-1],[2:n n])+u([2:m m],[1 1:n-1]);
	uxy = (Dp-Dm)/4;
    % compute flow
    %eq2 is introduced to make the solution unconditionally stable
    %curv is calculate according to
    %https://ww2.mathworks.cn/en/company/newsletters/articles/applying-modern-pde-techniques-to-digital-image-restoration.html
    %part "Factoring in Gradients: Fourth-Order Total Variation Equation"
    Num = uxx.*(ep2+uy.^2)-2*ux.*uy.*uxy+uyy.*(ep2+ux.^2);
    Den = (ep2+ux.^2+uy.^2).^(3/2);
    curv = Num./Den;
    lapcurvbar = Denominator.*fft2(curv);
    
    %it's the Cahn-Hilliard equation expressed in fourier space
    %this time step scheme comes from
    %file:///E:/%E6%B5%8F%E8%A7%88%E5%99%A8/Cahn-Hilliard_Inpainting_and_a_Generalization_for_.pdf
    %page 19
    ubar = ((1+c2*dt+dt*c1*Denominator.^2).*ubar - dt*lapcurvbar ...
        + dt*(lu0bar-lubar))./(1+dt*c2+dt*c1*Denominator.^2);
    
    %ifft is the reverse of fourier transform
    u = real(ifft2(ubar));
    it=it+1;
    figure(2);
    imagesc(u,clims); axis image; axis off; colormap(gray);
    title(['Inpainted image as solution to the fourth0-order total '...
        'variation flow after ' num2str(it) ' timesteps']);
    pause(0.01)
    lubar=fft2(lambda.*u);
end
