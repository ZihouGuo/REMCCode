function [UU,V,A,Z,iter,obj] = missingalgo_qp(X,Y,lambda,numanchor,ind)
%[F,V,A,Z,iter,obj] = missingalgo_qp(X1,gt,lambda1,lambda2,ind); % X,Y,lambda,d,numanchor
% m      : the number of anchor. the size of Z is m*n.
% lambda : the hyper-parameter of regularization term.

% X      : n*di

%% initialize
maxIter = 50 ; % the number of iterations

m = numanchor;
numclass = length(unique(Y));
numview = length(X);
numsample = size(Y,1);
% Z = zeros(m,numsample); 
% Z(:,1:m) = eye(m);
Z = zeros(m,numsample,numview); 
for iiv = 1:numview
Z(:,1:m,iiv) = eye(m);
end
Y = zeros(m,numsample,numview); 
Phi = zeros(m,numsample,numview); % ÕÅÁ¿
mu = numsample*1e-6;  
max_mu = 1e10;
missingindex = constructA(ind);
for i = 1:numview
    di = size(X{i},1); 
    A{i} = zeros(di,m); %Ãªµã¾ØÕó
end

alpha = ones(1,numview)/numview;

flag = 1;
iter = 0;
%%
while flag
    iter = iter + 1;

    %% optimize Ai
    for ia = 1:numview
       %part1 = X{ia} * Z';
        part1 = X{ia} * Y(:,:,ia)';
        [Unew,~,Vnew] = svd(part1,'econ');
        A{ia} = Unew*Vnew';
    end
    rho=1.2;
    sX= [m,numsample,numview];
    %% optimize Z
    temp          = Y - Phi/mu;
    [vecz, ~]     = wshrinkObj(temp(:),1/rho,sX,0,3);
    Z             =  reshape(vecz, sX);
    %% optimize Y
%         C1 = zeros(:,:,numview);
%         C2 = zeros(:,:,numview);
        C1 =0;
        C2 =0;
        for a=1:numview
            C1=lambda*alpha(a)^2*ind(:,a)'+0.5*mu * ones(1,numsample); 
            C2=lambda*alpha(a)^2 * A{a}'*X{a}+0.5*mu*Z(:,:,a)+0.5*Phi(:,:,a);
            
            for ii=1:numsample
               idx = 1:numanchor;
               ut = C2(idx,ii)./C1(ii);
               %Z(idx,ii) = EProjSimplex_new(ut');
               Y(idx,ii,a)=EProjSimplex_new(ut');
           end
        end

    %% optimize alpha
    M = zeros(numview,1);
    for iv = 1:numview
%         M(iv) = norm( X{iv} - A{iv} * (Z.*repmat(missingindex{iv},m,1)),'fro')^2;
          M(iv) = norm( X{iv} - A{iv} * (Y(:,:,iv).*repmat(missingindex{iv},m,1)),'fro')^2;
    end
    Mfra = M.^-1;
    Q = 1/sum(Mfra);
    alpha = Q*Mfra;

   %%
    Phi         = Phi + mu*(Z - Y);
    chag2       = max(abs(Z(:) - Y(:)));
    
    %%
    term1 = 0;
    term2 = 0;
    term3=0;
    for iv = 1:numview
        term1 = term1 + alpha(iv)^2 * norm(X{iv} - A{iv} * (Z(:,:,iv).*repmat(missingindex{iv},m,1)),'fro')^2;
        term3=term3+Z(:,:,iv);
    end
    term3=term3/numview;
    term2 = lambda * norm(term3,'fro')^2;
    obj(iter) = term1+term2;
    
    
    if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-3 || iter>maxIter || obj(iter) < 1e-10)
%         [UU,~,V]=svd(Z','econ');
        [UU,~,V]=svd(term3','econ');
        UU = UU(:,1:numclass);
        flag = 0;
    end
end
         
         
    
