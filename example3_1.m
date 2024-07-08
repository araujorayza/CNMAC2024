clear all;
clc;
close all;

%system
A{1}=[-10 30;0 -20];
A{2}=[-10 -30;0 -20];
n=2;

h{1} = @(x1,x2) (x1.^2+x2.^2)/50;
h{2} = @(x1,x2) 1-(x1.^2+x2.^2)/50;

dh{1} = @(x1,x2) [ 2.*x1/50,  2.*x2/50];
dh{2} = @(x1,x2) [-2.*x1/50, -2.*x2/50];

G=[1];


Rset=1:n;

Z=-5:0.05:5;

%LMI calculations
LMIS=[];
for j=G
    P{j} = sdpvar(n,n,'symmetric');
    R{j} = sdpvar(n,n,'full');
    L{j} = sdpvar(n,n,'full');
end

% b=0.35;
% l=0.2;

lambda=.8;
l=0.3
% l=0.15
% sdpvar l

for j=Rset
    for k=G
        Upsilon{k,j} = [L{k}*A{j}+A{j}'*L{k}'+lambda*P{k},   (P{k}-L{k}'+R{k}*A{j})', zeros(n,1);
                        P{k}-L{k}'+R{k}*A{j},        -R{k}-R{k}', zeros(n,1);
                        zeros(1,n), zeros(1,n), -lambda*l];
        LMIS = [LMIS, Upsilon{k,j} <= 0];
    end
end

% Less conservative LMIs
for k=G
    LMIS = [LMIS, Upsilon{k,k} <= 0];
end

for j=G
    for k=G
        if(k~=j)
            LMIS = [LMIS, Upsilon{k,j}+Upsilon{j,k} <= 0];
        end
    end
end

for j=setdiff(Rset,G)
    for k=G
        if(k~=j)
            LMIS = [LMIS, Upsilon{k,j} <= 0];
        end
    end
end

LMIS = [LMIS, lambda >= 0];

opts=sdpsettings;
opts.solver='sedumi';
opts.verbose=0;

sol = solvesdp(LMIS,[],opts);
p=min(checkset(LMIS));
if p > 0
    for k = G
        P{k} = double(P{k})
        R{k} = double(R{k})
        L{k} = double(L{k})
        lambda = double(lambda)
        l=double(l)
    end
else
    display('Infeasible')
    P=[];
end

% Set estimation
V = @(x1,x2) sum(arrayfun(@(k) [x1;x2]'*h{k}(x1,x2)*P{k}*[x1;x2],G));
hdot = @(x1,x2,k) sum(arrayfun(@(j) dh{k}(x1,x2)*h{j}(x1,x2)*A{j}*[x1;x2],Rset));
Dset = @(x1,x2) sum(arrayfun(@(k) [x1;x2]'*hdot(x1,x2,k)*P{k}*[x1;x2],G));

%calculate V and D
meshPoints=500;
tol=10/meshPoints;
x = linspace(-5,5,meshPoints);
y = linspace(-5,5,meshPoints);
[X,Y]=meshgrid(x,y);
for i=1:length(x)
    for j = 1:length(y)
        Ve(i,j) = V(X(i,j),Y(i,j));
        De(i,j) = Dset(X(i,j),Y(i,j));
    end
end

%calculate b
b=min([min(Ve(:,1)), min(Ve(:,end)), min(Ve(1,:)), min(Ve(end,:))])
b=fix(b*1e2)/1e2;

figure(1);
[~,c]=contour(X,Y,Ve,linspace(0,b,5),'r','ShowText','on','DisplayName','V')
hold on
[~,d]=contour(X,Y,De,[0,fix(max(max(De))*1e2)/1e2],'b','ShowText','on','DisplayName','D')
legend;
%from the graph
%l = 0.175;



