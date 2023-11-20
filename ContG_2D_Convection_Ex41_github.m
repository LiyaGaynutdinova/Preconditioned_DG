function  Fce_hom_1
% Continuous galerkin, diffusion-conduction-reaction, piecewise linear functions
% (0,L1) x (0,L2)

% close all
clear all;
Uniform = 1; % 1=uniform or 0=nonuniform (so far 1 is working only)
toler = 1e-8;
max_steps = 5000;
% sigma = 10;
L1 = 1; L2 = 1; 
plot_fig = 1; % do we wish to plot graphs? 0 or 1

Nintervals = 20;

NNN = ones(2,1)*Nintervals; 
N1 = NNN(1); N2 = NNN(2); % numbers of intervals
x1 = linspace(0,L1,N1+1);
x2 = linspace(0,L2,N2+1);
if (Uniform==1) 
    h1 = L1/N1; h2 = L2/N2; 
    Ninner = (N1-1)*(N2-1);
end;
[X1,X2] = meshgrid(x1,x2);
[nodes,bnodes] = Nodes(N1,N2,L1,L2); % all nodes, including boundary nodes
elem = Elem(N1,N2,nodes);
Nele = max(size(elem))
elemNew = delaunay(X1,X2);
NeleNew = max(size(elemNew))
for k = 1:NeleNew
    pom = elemNew(k,:);
    pom = sort(pom);
    elemNew(k,:) = pom; 
end;
Nnod = (N1+1)*(N2+1)

% continuous Galerkin:
A = zeros(Nnod); % stiffness matrix + sym Nav
Nav = A; % convection matrix - non-sym;
AP = A;
B = zeros(Nnod,1); % rhs
D1 = zeros(Nele,Nnod); % 2 matrices of derivatives in quadrature points (dx and dy)
D2 = zeros(Nele,Nnod);
upp1 = zeros(Nele,1);
low1 = zeros(Nele,1);
uppNod = zeros(Nnod,1)-1;
lowNod = zeros(Nnod,1)+1e10;
uppNodIm = zeros(Nnod,1)-1;
imag_max = zeros(Nele,1);
APcoord = zeros(Nele,3);
for kk = 1:Nele
    Apom = zeros(3);
    APpom = zeros(3);
    Navpom = zeros(3);
    kde = elem(kk,1:3);
    nod = nodes(kde,:);
    dd = abs(det([nod(1,:)-nod(3,:);nod(2,:)-nod(3,:)]))/2;
    nod1 = [nod,ones(3,1)];
    der = inv(nod1);
    der = der(1:2,:); % d/dx, d/dy for 3 basis functions
    g = sum(nod)/3; % mass center
    aa = FA(g(1),g(2)); 
    aap = FAP(g(1),g(2));
    ff = FF(g(1),g(2));  
    rr = FR(g(1),g(2));
    rrp = FRP(g(1),g(2));    
    pom8 = [kk,kk+1*Nele];
    % derivatives [dxu1,0,dyu1]' and [0,dyu2,dxu2]' in all quad.points:    
    derr = [der(1,:);der(2,:)]; 
    pom = derr'*aa*derr;
    pomp = derr'*aap*derr;
    Apom = pom*dd + rr*[2,1,1;1,2,1;1,1,2]*h1*h2/24;  % ? correct
    APpom = pomp*dd + rrp*[2,1,1;1,2,1;1,1,2]*h1*h2/24;  % ? correct
    A(kde,kde) = A(kde,kde) + Apom;
    AP(kde,kde) = AP(kde,kde) + APpom;
    B(kde) = B(kde)+ff*dd/6;
    
    bb = Fb(g(1),g(2)) ;  
    pomNav = derr'*bb*[1,1,1]*h1*h2/2;   % convection term 

    pomNav1 = (pomNav+pomNav')/2; 
    pomNav2 = (pomNav-pomNav')/2;
    
    [u1,c2] = eig(pomNav1);
    pomNavSpec = pomNav1;

    Navpom = pomNav2;
    Apom = Apom + pomNav1*0;
    APpom = APpom + pomNavSpec*0;
    A(kde,kde) = A(kde,kde) + pomNav1*0;
    AP(kde,kde) = AP(kde,kde) + pomNavSpec*0;
    Nav(kde,kde) = Nav(kde,kde) + pomNav2;
       
    pom8 = sort(eig(Apom));  
    pom8p = sort(eig(APpom));  
%     if (min(pom8)<-0.0001  || min(pom8p)<-0.0001) 
    if (min(pom8p)<-0.0001)     
       disp('=========== chyba ============');  
       [pom8,pom8p]
       APpom
       return; end;   
    pom9 = pinv(APpom)*Apom;   
    p9 = sort(real(eig(pom9)));
    pom10 = max(imag(eig(pinv(APpom)*Navpom)));
     if (abs(p9(1))<0.00001) low1(kk) = p9(2);
    else low1(kk) = p9(1); end;
%     low1(kk) = p9(1);
    upp1(kk) = p9(end);
    for jj3 = 1:3
    lowNod(kde(jj3)) = min([p9(1),lowNod(kde(jj3))]);
    uppNod(kde(jj3)) = max([p9(end),uppNod(kde(jj3))]);
    uppNodIm(kde(jj3)) = max([pom10,uppNodIm(kde(jj3))]);
    end;
    imag_max(kk) = pom10;
    APcoord(kk,:) = kde;
end;
low2 = zeros(Nnod,1)+1e8;
upp2 = zeros(Nnod,1)-100;
imag_max2 = zeros(Nnod,1)-100;
APc1 = zeros(Nnod,3);
APc2 = APc1;
for kk = 1:Nele % 
    for j = 1:3
        s = elem(kk,j);
        if (low1(kk)<low2(s)) low2(s)=low1(kk); APc1(s,:) = APcoord(kk,:); end;       
        if (upp1(kk)>upp2(s)) upp2(s)=upp1(kk); APc2(s,:) = APcoord(kk,:); end;
        if (imag_max(kk)>imag_max2(s)) imag_max2(s)=imag_max(kk); end;
    end;        
end;
for k = Nnod:-1:1 
if (bnodes(k)==1) A(k,:) = []; A(:,k) = []; B(k) = []; 
    AP(k,:) = []; AP(:,k) = []; 
    upp2(k) = []; low2(k) = []; imag_max2(k) = [];    
    APc1(k,:) = [];  APc2(k,:) = [];
    Nav(k,:) = []; Nav(:,k) = [];   
end;
end;
% better lower and upper bounds: 
 y = ContGalBounds(elem,Nele,Nnod,bnodes,low1,upp1,Ninner);
 low3 = y(1:Ninner,1);
 upp3 = y(1:Ninner,2);

%solution U:
U = (A+Nav)\B;
U1 = reshape(U,N1-1,N2-1);
U2 = zeros(N1+1,N2+1);
U2(2:N1,2:N2) = U1';
if (plot_fig==1)
    subplot(2,2,1); % solution
    surf(X1,X2,U2);
    view(15,15);
    title('solution');
end;

eigA = sort(real(eig(A)));
eigPA = sort(real(eig(inv(AP)*A)));
conditionA = max(eigA)/min(eigA)
conditionPA = max(eigPA)/min(eigPA)
conditionPAbound = max(upp2)/min(low2)

if (plot_fig==1)
    subplot(2,2,2); % real part of the spectrum of P^(-1)A
    cla; hold on;
    plot(sort(upp2),'b');
    plot(eigPA,'k.');
    plot(sort(low2),'r');
    plot(low3(1:(N1-1)*(N2-1)),'k--'); % alternative bound
    plot(upp3((N1-1)*(N2-1):-1:1),'k--'); % alternative bound
    title('real parts of eigenvalues');
    
    subplot(2,2,3);
    cla; hold on;
    spy(AP);
    title('preconditioner');
end;

p = eig(A+Nav);
p2 = eig(inv(AP)*(A+Nav));
my1 = max(imag_max2);
mx1 = min(low2);
mx2 = max(upp2);
eigAP = sort(real(eig(inv(AP)*A)));
eigBP = sort(imag(eig(inv(AP)*Nav)));
my1 = max(eigBP);
mx1 = min(eigAP);
mx2 = max(eigAP);

if (plot_fig==1)
    subplot(2,2,4);
    cla; hold on;
    plot(real(p),imag(p),'k.');
    plot(real(p2),imag(p2),'r.');
    title('eigenvales and bounds');
    plot([mx1,mx2,mx2,mx1,mx1],[-my1,-my1,my1,my1,-my1],'k--'); % worse but feasible
    plot([mx1,mx2,mx2,mx1,mx1],[-my1,-my1,my1,my1,-my1],'b'); % better but unfeasible
end;

maxBeig = max(imag(eig(Nav))) % max imag part of the eigenvalues of B
maxPBeig = my1 % max imag part of the eigenvalues of P(-1)B
maxPBeig_estim = max(imag_max2) % bound to imag parts of the eigenvalues of P(-1)B

n = max(size(A));
tic;
[x,flagA,relres,iterA] = gmres(A+Nav,B,[],toler,n);
tA = toc;
tic;
[x,flagPA2,relres,iterPA2] = gmres(A+Nav,B,[],toler,n,AP);
tPA2 = toc;
flags12 = [flagA,flagPA2]
iter12 = [iterA,iterPA2]
times12 = [tA,tPA2];

 

%============================================================
%============================================================
%============================================================
%============================================================

function y = FA(x1,x2)
y = [20-2*x2,0;0,3-2*x1];

function y = FAP(x1,x2)
y = [19,0;0,2];
% y = [1,0;0,1];

function y = FR(x1,x2) % reaction coefficient
y = 10;

function y = FRP(x1,x2) % reaction coefficient
y = 10;

function y = Fb(x1,x2) % convection coefficient
y = [-x2;x1]*10;

function y = FF(x1,x2)
y = 10;

function [nodes,bnodes] = Nodes(N1,N2,L1,L2)
NN = (N1+1)*(N2+1);
h1 = L1/N1; h2 = L2/N2; 
nodes = zeros(NN,2);
bnodes = zeros(NN,1);
j = 1;
for k2 = 1:N2+1
    for k1 = 1:N1+1        
            nodes(j,:) = [(k1-1)*h1,(k2-1)*h2];
            if (k1==1 || k1==N1+1 || k2==1 || k2==N2+1)
            bnodes(j) = 1; end;
            j = j+1;
    end;
end;

function elem = Elem(N1,N2,nodes)
NN = N1*N2;
elem = zeros(2*NN,3);
j = 1;
for k2 = 1:N2
    for k1 = 1:N1
        p = (k2-1)*(N1+1)+k1; % lower left corner
            elem(j,:) = [p,p+1,p+N1+2]; 
            elem(j+1,:) = [p,p+N1+1,p+N1+2];
            j = j+2;
    end;
end;

function y = ContGalBounds(elem,Nele,Nnod,bnodes,low1,upp1,Ninner)
% better lower bounds:
low3 = zeros(Nnod,1);
elem_eig = elem;
for kk = 1:Nele
    for jj = 1:3
        if (bnodes(elem_eig(kk,jj))>0) elem_eig(kk,jj)=0; end;
    end;
    if (sum(elem_eig(kk,:))==0) low1(kk)=100000; upp1(kk)=-1; end;
end;
for jj = 1:Nnod       
    [vL,eL] = min(low1);
    low3(jj) = vL;
    pom = elem_eig(eL,:);
    [n1,c] = max(pom);
    for kk = 1:Nele
        for mm = 1:3
            if (elem_eig(kk,mm)==n1) elem_eig(kk,mm)=0; end            
        end;
        if (sum(elem_eig(kk,:))==0) low1(kk)=1e9; end;
    end;
%     [low1,elem_eig];
end;
% better lower bounds:
upp3 = zeros(Nnod,1);
elem_eig = elem;
for kk = 1:Nele
    for jj = 1:3
        if (bnodes(elem_eig(kk,jj))>0) elem_eig(kk,jj)=0; end;
    end;
    if (sum(elem_eig(kk,:))==0) upp1(kk)=-10000000; upp1(kk)=-1; end;
end;
ccc = [upp1,elem_eig];
for jj = 1:Nnod       
    [vU,eU] = max(upp1);
    upp3(jj) = vU;
    pom = elem_eig(eU,:);
    [n1,c] = max(pom);
    for kk = 1:Nele
        for mm = 1:3
            if (elem_eig(kk,mm)==n1) elem_eig(kk,mm)=0; end            
        end;
        if (sum(elem_eig(kk,:))==0) upp1(kk)=-1e9; end;
    end;
%     [upp1,elem_eig];
end;
y = [low3,upp3];

