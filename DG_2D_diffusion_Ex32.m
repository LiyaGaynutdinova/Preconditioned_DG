function  Fce_abc
% Discontinuous galerkin, diffusion-reaction, piecewise linear functions
% (comparison with continuous Lagrange FEM)
% (0,L1) x (0,L2)
% boundary edges icluded !

% close all
clear all;
clc;

Uniform = 1; % 1=uniform or 0=nonuniform (so far 1 is working only)
toler = 1e-6;
max_steps = 50000; % for CG
L1 = 1; L2 = 1; 
c_sigma = 2;  % sigma_p for the original problem ( >=1)
c_sigma_p = 2; % sigma_p for preconditioning problem ( >=1)

NNNNN = 10 %%% = 1/h

NNN = ones(2,1)*NNNNN; % number of intervals in every direction (2D)

N1 = NNN(1); N2 = NNN(2); % number of intervals
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
Nnod = (N1+1)*(N2+1)
Nvox = N1*N2;  % number of pixels
more2few = More2Few(elem,Nele); % structure: elem, node_few
edge2elem = Edge2Elem(elem,Nele); 
% structure: n1,n2,e1,e2,n11,n12,n13,n21,n22,n23,m11,m12,m13,m21,m22,m23
pom = size(edge2elem);
Nedg = pom(1);
Ndof = 3*Nele;

% continuous Galerkin: ======================================
A = zeros(Nnod); % stiffness matrix
AP = A;
B = zeros(Nnod,1); % rhs
D1 = zeros(Nele,Nnod); % 2 matrices of derivatives in quadrature points (dx and dy)
D2 = zeros(Nele,Nnod);
upp1 = zeros(Nele,1);
low1 = zeros(Nele,1);
APcoord = zeros(Nele,3);
for kk = 1:Nele
    Apom = zeros(3);
    APpom = zeros(3);
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
    Apom = pom*dd + rr*[2,1,1;1,2,1;1,1,2]*h1*h2/24;  % ?correct
    APpom = pomp*dd + rrp*[2,1,1;1,2,1;1,1,2]*h1*h2/24;  % ?correct
    A(kde,kde) = A(kde,kde) + Apom;
    AP(kde,kde) = AP(kde,kde) + APpom;
    B(kde) = B(kde)+ff*dd/6;
    pom8 = sort(eig(Apom));  
    pom8p = sort(eig(APpom));  
    if (min(pom8)<-0.0001  || min(pom8p)<-0.0001) 
       disp('=========== chyba ============');  
       [pom8,pom8p]
       return; end;   
    pom9 = pinv(APpom)*Apom;   
    p9 = sort(real(eig(pom9)));
    if (abs(p9(1))<0.00001) low1(kk) = p9(2);
    else low1(kk) = p9(1); end;
    upp1(kk) = p9(end);
    APcoord(kk,:) = kde;
end;
low2 = zeros(Nnod,1)+1e8;;
upp2 = zeros(Nnod,1)-100;
APc1 = zeros(Nnod,3);
APc2 = APc1;
for kk = 1:Nele % 
    for j = 1:3
        s = elem(kk,j);
        if (low1(kk)<low2(s)) low2(s)=low1(kk); APc1(s,:) = APcoord(kk,:); end;       
        if (upp1(kk)>upp2(s)) upp2(s)=upp1(kk); APc2(s,:) = APcoord(kk,:); end;
    end;        
end;
for k = Nnod:-1:1 
if (bnodes(k)==1) A(k,:) = []; A(:,k) = []; B(k) = []; 
    AP(k,:) = []; AP(:,k) = []; 
    upp2(k) = []; low2(k) = [];
    APc1(k,:) = [];  APc2(k,:) = [];
end;
end;
% possibly better lower and upper bounds: 
 y = ContGalBounds(elem,Nele,Nnod,bnodes,low1,upp1,Ninner);
 low3 = y(1:Ninner,1);
 upp3 = y(1:Ninner,2);

%solution U:
U = A\B;
U1 = reshape(U,N1-1,N2-1);
U2 = zeros(N1+1,N2+1);
U2(2:N1,2:N2) = U1';
% subplot(2,2,1);
% surf(X1,X2,U2); 
% title('Cont G solution')
% view(15,15);
% subplot(2,2,2);
% cla; hold on;
eigA = sort(real(eig(A)));
eigPA = sort(real(eig(inv(AP)*A)));
conditionA = max(eigA)/min(eigA)
conditionPA = max(eigPA)/min(eigPA)
% plot(sort(upp2),'b');
% plot(eigPA,'k.');
% plot(sort(low2),'r');
% plot(low3(1:(N1-1)*(N2-1)),'k--');
% plot(upp3((N1-1)*(N2-1):-1:1),'k--');
% title('Cont G: eig of P^{-1}A and bounds');
conditionPAbound = max(upp2)/min(low2)

[x,flag,relres,iterPG]  = cgs(A,B,toler,max_steps,AP);
[x,flag,relres,iterG]  = cgs(A,B,toler,max_steps);
iterationsGandPG = [iterG,iterPG]



% discontinuous Galerkin: =============================================
ADG = zeros(Ndof); % DOFs - nodal values of disc. lin. function
APDG = ADG;
BDG = zeros(Ndof,1);
lower = zeros(Nedg,1);
upper = zeros(Nedg,1);
minAeig = 1e10;
maxAeig = -1;
for kk = 1:Nedg    
    Apom = zeros(6);
    APpom = zeros(6);
    p = edge2elem(kk,:);

% structure of a row: n1,n2,e1,e2,n11,n12,n13,n21,n22,n23,m11,m12,m13,m21,m22,m23
% (2 nodes, 2 elements (left, right), 
% 3 unique global numbers of nodes (left), 3 unique global numbers of nodes (right), 
% 3 longer global numbers of nodes (left), 3 longer global numbers of nodes (right))

    if (p(4)==0) % 1/3 of diffusion, reaction and F and both boundary integrals
        e1 = p(3);
        dof1 = p([11,12,13]); % nodesDG in first elem
        n1 = p([5,6,7]); % nodes in first elem
        point1 = nodes(n1,:);
        cen1 = sum(point1)/3;
        % element1 integral 1/3: 
        dd1 = abs(det([point1(1,:)-point1(3,:);point1(2,:)-point1(3,:)]))/2;
        nod1 = [point1,ones(3,1)];
        der1 = inv(nod1);
        derr1 = der1(1:2,:); % d/dx, d/dy for 3 basis functions
        aa1 = FA(cen1(1),cen1(2)); 
        ff1 = FF(cen1(1),cen1(2));
        pp1 = FAP(cen1(1),cen1(2));         
        pom1 = derr1'*aa1*derr1;  
        ADG(dof1,dof1) = ADG(dof1,dof1) + pom1*dd1/3;
        APDG(dof1,dof1) = APDG(dof1,dof1) + derr1'*pp1*derr1*dd1/3;
        BDG(dof1) = BDG(dof1) + ff1*dd1/6/3;
         Apom = pom1*dd1/3;
         APpom = derr1'*pp1*derr1*dd1/3;
        pom9 = pinv(APpom)*Apom;
        p9 = sort(real(eig(pom9)));
        lower(kk) = p9(2);
        upper(kk) = p9(end);
        
        % mass matrices: [2,1,1;1,2,1;1,1,2]*h^2/24;  % only ONE THIRD!
        rr1 = FR(cen1(1),cen1(2));    
        prr1 = FRP(cen1(1),cen1(2));    
        ADG(dof1,dof1) = ADG(dof1,dof1) + rr1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
        APDG(dof1,dof1) = APDG(dof1,dof1) + prr1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;    
        Apom(1:3,1:3) = Apom(1:3,1:3) + rr1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;    
        APpom(1:3,1:3) = APpom(1:3,1:3) + prr1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
         
        % edge integral:
        p1 = nodes(n1(1),:);
        p2 = nodes(n1(2),:);    
        de = norm(p2-p1); % edge length
        nor1 = p2 - p1;
        nor1 = [nor1(2),-nor1(1)];    
        nor1 = nor1/norm(nor1); % unit normal (directing to right)    
        flux1 = derr1'*aa1*nor1';       
        flux1p = derr1'*pp1*nor1';
        pom1 = [flux1]/2;    
        pom1p = [flux1p]/2;
        pom2 = [1,1,0]/2;
        flux = -(pom1*pom2+pom2'*pom1')*de; % co zde, + nebo - ????
        fluxP = -(pom1p*pom2+pom2'*pom1p')*de;  
        ADG([dof1],[dof1]) = ADG([dof1],[dof1]) + flux;
        APDG([dof1],[dof1]) = APDG([dof1],[dof1]) + fluxP; % !!!
        Apom = Apom + flux;    
        APpom = APpom + fluxP;  % !!!
 
        % penalty term:
        lam_aL = real(eig(aa1));        
        lam_pL = real(eig(pp1));        
        sigmaloc = 12/de*(max(lam_aL)^2/min(lam_aL) ) *  c_sigma;
        sigmaloc_p = 12/de*(max(lam_pL)^2/min(lam_pL) )  * c_sigma_p;    

        pom7 = [2,1;1,2]*de/6; % integrals of [[u]][[v]] 

        APDG([dof1(1:2)],[dof1(1:2)]) = APDG([dof1(1:2)],[dof1(1:2)])+pom7*sigmaloc_p; % !!!
        Apom([1,2],[1,2]) =  Apom([1,2],[1,2]) + pom7*sigmaloc;
        APpom([1,2],[1,2]) =  APpom([1,2],[1,2]) + pom7*sigmaloc_p;% !!!   
        
        % check of positivity:
        pom8 = sort(real(eig(Apom)));  
        pom8p = sort(real(eig(APpom)));  
        if (min(pom8)<-0.0000001  || min(pom8p)<-0.000001) 
            disp('=========== error ============'); % wrong choice of c_sigma
            [pom8,pom8p]
            kk
        return; end;
        if (pom8(1)<minAeig) minAeig = pom8(1); end;
        if (pom8(end)>maxAeig) maxAeig = pom8(end); end;
        % eigen-bounds preparing: 
        pom9 = pinv(APpom)*Apom;
        p9 = sort(real(eig(pom9)));
        if (abs(p9(1))<0.000001) lower(kk) = p9(2);
        else lower(kk) = p9(1); end;
%       lower(kk) = p9(1);
        upper(kk) = p9(end);   
        continue;   
    end;
    n1 = edge2elem(kk,[5,6,7]); %elem(e1,:);   
    n2 = edge2elem(kk,[8,9,10]); %elem(e2,:);
    dof1 = edge2elem(kk,11:13); % 
    dof2 = edge2elem(kk,14:16); % 
    sdof1 = sort(dof1);
    sdof2 = sort(dof2);
    point1 = nodes(n1,:);
    point2 = nodes(n2,:);
    cen1 = sum(point1)/3;
    cen2 = sum(point2)/3;     
    % element1 integral 1/3:    
    dd1 = abs(det([point1(1,:)-point1(3,:);point1(2,:)-point1(3,:)]))/2;
    nod1 = [point1,ones(3,1)];
    der1 = inv(nod1);
    derr1 = der1(1:2,:); % d/dx, d/dy for 3 basis functions
    aa1 = FA(cen1(1),cen1(2)); 
    pp1 = FAP(cen1(1),cen1(2));
    ff1 = FF(cen1(1),cen1(2)); 
    % derr1 = [der1(1,:);der1(2,:)]; 
    pom1 = derr1'*aa1*derr1;
    pom1p = derr1'*pp1*derr1;
    ADG(dof1,dof1) = ADG(dof1,dof1) + pom1*dd1/3;
    APDG(dof1,dof1) = APDG(dof1,dof1) + pom1p*dd1/3;
    BDG(dof1) = BDG(dof1) + ff1*dd1/6/3;   
    % element2 integral 1/3:    
    dd2 = abs(det([point2(1,:)-point2(3,:);point2(2,:)-point2(3,:)]))/2;
    nod2 = [point2,ones(3,1)];
    der2 = inv(nod2);
    derr2 = der2(1:2,:); % d/dx, d/dy for 3 basis functions
    aa2 = FA(cen2(1),cen2(2)); 
    pp2 = FAP(cen2(1),cen2(2)); 
    ff2 = FF(cen2(1),cen2(2));
    % derr2 = [der2(1,:);der2(2,:)]; 
    pom2 = derr2'*aa2*derr2;
    pom2p = derr2'*pp2*derr2;
    ADG(dof2,dof2) = ADG(dof2,dof2) + pom2*dd2/3;
    APDG(dof2,dof2) = APDG(dof2,dof2) + pom2p*dd2/3;
    BDG(dof2) = BDG(dof2) + ff2*dd2/6/3;    
    Apom(1:3,1:3) = pom1*dd1/3;
    Apom(4:6,4:6) = pom2*dd2/3;
    APpom(1:3,1:3) = pom1p*dd1/3;
    APpom(4:6,4:6) = pom2p*dd2/3;  
   
  % mass matricces: [2,1,1;1,2,1;1,1,2]*h^2/24;  % only ONE THIRD!
    rr1 = FR(cen1(1),cen1(2));
    rr2 = FR(cen2(1),cen2(2));
    prr1 = FRP(cen1(1),cen1(2));
    prr2 = FRP(cen2(1),cen2(2));
    ADG(dof1,dof1) = ADG(dof1,dof1) + rr1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
    ADG(dof2,dof2) = ADG(dof2,dof2) + rr2*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
    APDG(dof1,dof1) = APDG(dof1,dof1) + prr1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
    APDG(dof2,dof2) = APDG(dof2,dof2) + prr2*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
    Apom(1:3,1:3) = Apom(1:3,1:3) + rr1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
    Apom(4:6,4:6) = Apom(4:6,4:6) + rr2*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
    APpom(1:3,1:3) = APpom(1:3,1:3) + prr1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;
    APpom(4:6,4:6) = APpom(4:6,4:6) + prr2*[2,1,1;1,2,1;1,1,2]*h1*h2/24/3;    
       
  % edge integral:
    p1 = nodes(n1(1),:);
    p2 = nodes(n1(2),:);    
    de = norm(p2-p1); % edge length
    nor1 = p2 - p1;
    nor1 = [nor1(2),-nor1(1)];    
    nor1 = nor1/norm(nor1); % unit normal (directing to right)
    
    flux1 = derr1'*aa1*nor1';
    flux2 = derr2'*aa2*nor1';
    flux1p = derr1'*pp1*nor1';
    flux2p = derr2'*pp2*nor1';
    pom1 = [flux1;flux2]/2;    
    pom1p = [flux1p;flux2p]/2;
    pom2 = [1,1,0,-1,-1,0]/2;
    flux = -(pom1*pom2+pom2'*pom1')*de; % co zde, + nebo - ????
    fluxP = -(pom1p*pom2+pom2'*pom1p')*de;    
    
    ADG([dof1,dof2],[dof1,dof2]) = ADG([dof1,dof2],[dof1,dof2]) + flux;
    APDG([dof1,dof2],[dof1,dof2]) = APDG([dof1,dof2],[dof1,dof2]) + fluxP; % !!!
    Apom = Apom + flux;    
    APpom = APpom + fluxP;  

  % penalty term:
  lam_aL = real(eig(aa1));
  lam_aR = real(eig(aa2));
  lam_pL = real(eig(pp1));
  lam_pR = real(eig(pp2));
  sigmaloc = 6/de*(max(lam_aL)^2/min(lam_aL) + max(lam_aR)^2/min(lam_aR)) * c_sigma;
  sigmaloc_p = 6/de*(max(lam_pL)^2/min(lam_pL) + max(lam_pR)^2/min(lam_pR))  * c_sigma_p;

    pom7 = [2,1,-2,-1;1,2,-1,-2;-2,-1,2,1;-1,-2,1,2]*de/6; % integrals of [[u]][[v]]   

    ADG([dof1(1:2),dof2(1:2)],[dof1(1:2),dof2(1:2)]) = ...
        ADG([dof1(1:2),dof2(1:2)],[dof1(1:2),dof2(1:2)])+pom7*sigmaloc;
    APDG([dof1(1:2),dof2(1:2)],[dof1(1:2),dof2(1:2)]) = ...
        APDG([dof1(1:2),dof2(1:2)],[dof1(1:2),dof2(1:2)])+pom7*sigmaloc_p; % !!!
    Apom([1,2,4,5],[1,2,4,5]) =  Apom([1,2,4,5],[1,2,4,5]) + pom7*sigmaloc;
    APpom([1,2,4,5],[1,2,4,5]) =  APpom([1,2,4,5],[1,2,4,5]) + pom7*sigmaloc_p;% !!!   
    
    % check of positivity:
    pom8 = sort(real(eig(Apom)));  
    pom8p = sort(real(eig(APpom)));  
    if (min(pom8)<-0.0000001  || min(pom8p)<-0.000001) 
       disp('=========== chyba ============'); 
       [pom8,pom8p]
       kk
       return; end;
   if (pom8(1)<minAeig) minAeig = pom8(1); end;
   if (pom8(end)>maxAeig) maxAeig = pom8(end); end;
   % eigen-bounds preparing: 
   pom9 = pinv(APpom)*Apom;
   p9 = sort(real(eig(pom9)));
   if (abs(p9(1))<0.00001) lower(kk) = p9(2);
   else lower(kk) = p9(1); end;
   if (abs(p9(2))<0.00001) lower(kk) = p9(3); end;   
%    lower(kk) = p9(1);
   upper(kk) = p9(end);   
end;

pom6 = zeros(Ndof,1);
for kk = 1:Nedg    
    if (bnodes(edge2elem(kk,5))==1) pom6(edge2elem(kk,11))=1; end;
    if (bnodes(edge2elem(kk,6))==1) pom6(edge2elem(kk,12))=1; end;
    if (bnodes(edge2elem(kk,7))==1) pom6(edge2elem(kk,13))=1; end;
   if (edge2elem(kk,4)>0)
    if (bnodes(edge2elem(kk,8))==1) pom6(edge2elem(kk,14))=1; end;
    if (bnodes(edge2elem(kk,9))==1) pom6(edge2elem(kk,15))=1; end;
    if (bnodes(edge2elem(kk,10))==1) pom6(edge2elem(kk,16))=1; end;
   end;     
end;
for kk = Ndof:-1:1
    if pom6(kk)==1
        ADG(:,kk) = [];
        ADG(kk,:) = [];
        APDG(:,kk) = [];
        APDG(kk,:) = [];
        BDG(kk) = [];
    end;
end;
% DG solution
UDG = ADG\BDG;
UDG1 = pom6*0;
UDG1(pom6<1) = UDG;

% subplot(2,2,3); % DG solution
% drawDG(UDG1,nodes,edge2elem) 
% view(15,15);
% title('DG solution');
PreDG = inv(APDG)*ADG;
eigADG = sort(real(eig(ADG)));
eigPDG = sort(real(eig(PreDG)));
conditionADG = eigADG(end)/eigADG(1)

PreDG = inv(APDG)*ADG;
conditionPADG = eigPDG(end)/eigPDG(1)

% subplot(2,2,4);
% cla; hold on;
% plot(eigPDG,'k.');
% title('DG: eig of P^{-1}A and bounds');

% eigen-bounds:
elower = zeros(Ndof,1)+1e8;
eupper = zeros(Ndof,1)-1;
for kk = 1:Nedg
    p = edge2elem(kk,:);
    c1 = lower(kk);
    c2 = upper(kk);    
    for jj = 11:13        
        if (c1<elower(p(jj))) elower(p(jj)) = c1; end;
        if (c2>eupper(p(jj))) eupper(p(jj)) = c2; end;        
    end;
    if (p(4)>0)
        for jj = 14:16        
        if (c1<elower(p(jj))) elower(p(jj)) = c1; end;
        if (c2>eupper(p(jj))) eupper(p(jj)) = c2; end;        
    end;
    end;
end;
for kk = Ndof:-1:1
    if pom6(kk)==1
        elower(kk) = [];
        eupper(kk) = [];
    end;
end;

% plot(sort(elower),'r');
% plot(sort(eupper),'b');
conditionPADGbound = max(eupper)/min(elower)
tic;
[x,flag,relres,iterPDG]  = cgs(ADG,BDG,toler,max_steps,APDG);
casPDG = toc;
tic;
[x,flag,relres,iterDG]  = cgs(ADG,BDG,toler,max_steps);
casDG = toc;
iterationsDGandPDG = [iterDG,iterPDG]
% timeDGandPDG = [casDG,casPDG]


%============================================================
%============================================================
%============================================================
%============================================================

function y = FA(x1,x2) % diffusion coefficient
y = [3.01+3*sin(x1*x2*pi)/1,0;0,1.01+sin(x1*x2*pi)/1];
% y = [1,0;0,1];
% if (x1>0.5) y = [5,0;0,10]; end;

function y = FAP(x1,x2) % preconditioning diffusion coefficient
y = [1,0;0,1];
% y = [3,0;0,1];

function y = FR(x1,x2) % reaction coefficient
y = 1;

function y = FRP(x1,x2) % reaction coefficient
y = 1;

function y = FF(x1,x2) % right hand side
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

function y = Edge2Elem(elem,Nele)
edge2elem = zeros(1,16); 
% structure: n1,n2,e1,e2,n11,n12,n13,n21,n22,n23,m11,m12,m13,m21,m22,m23
m = 1;
for kk = 1:Nele
    pp = elem(kk,:);
    p = [pp(1),pp(2)];
    moren = zeros(1,6);
    kam = 1;    
    for k3 = 1:m % already exists
        if ((edge2elem(k3,1)==p(1)&&(edge2elem(k3,2)==p(2)))) kam=2; break; end;
    end; 
    if (kam==1) 
        moren = [kk*3-2,kk*3-1,kk*3,0,0,0];
        newedge = [p(1),p(2),kk,0,pp(1),pp(2),pp(3),0,0,0,moren];
        edge2elem = [edge2elem;newedge]; 
        m = m+1;
    else edge2elem(k3,[4,8,9,10,14,15,16]) = [kk,pp(1),pp(2),pp(3),kk*3-2,kk*3-1,kk*3]; 
    end; 
    
    p = [pp(1),pp(3)];
    kam = 1;    
    for k3 = 1:m % already exists
        if ((edge2elem(k3,1)==p(1)&&(edge2elem(k3,2)==p(2)))) kam=2; break; end;
    end; 
    if (kam==1) 
        moren = [kk*3-2,kk*3,kk*3-1,0,0,0];
        newedge = [p(1),p(2),kk,0,pp(1),pp(3),pp(2),0,0,0,moren];
        edge2elem = [edge2elem;newedge]; 
        m = m+1;
    else edge2elem(k3,[4,8,9,10,14,15,16]) = [kk,pp(1),pp(3),pp(2),kk*3-2,kk*3,kk*3-1]; 
    end;   
    
    p = [pp(2),pp(3)];
    kam = 1;    
    for k3 = 1:m % already exists
        if ((edge2elem(k3,1)==p(1)&&(edge2elem(k3,2)==p(2)))) kam=2; break; end;
    end; 
    if (kam==1) 
        moren = [kk*3-1,kk*3,kk*3-2,0,0,0];
        newedge = [p(1),p(2),kk,0,pp(2),pp(3),pp(1),0,0,0,moren];
        edge2elem = [edge2elem;newedge]; 
        m = m+1;
    else edge2elem(k3,[4,8,9,10]) = [kk,pp(2),pp(3),pp(1),kk*3-1,kk*3,kk*3-2]; 
    end; 
end;
edge2elem(1,:) = [];
y = edge2elem;

function y = More2Few(elem,Nele)
more2few = zeros(3*Nele,2); % structure: elem, node_few
for kk = 1:Nele
    pom = elem(kk,:);
    w = [kk*3-2,kk*3-1,kk*3];
    more2few(w,:) = [kk,kk,kk;w]';    
end;
y = more2few;

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
% low1
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

function y = drawDG(UDG1,nodes,edge2elem) 
n = max(size(edge2elem));
cla; hold on;
for kk = 1:n
    n1 = edge2elem(kk,5:7);
    dof1 = edge2elem(kk,11:13);
    p1 = nodes(n1,:);
    x1 = p1(:,1);
    y1 = p1(:,2);    
    z1 = UDG1(dof1);
    surf([x1(1),x1(2);x1(3),x1(3)],[y1(1),y1(2);y1(3),y1(3)],[z1(1),z1(2);z1(3),z1(3)])
end;

