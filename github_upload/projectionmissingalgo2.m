function [UU,V,A,Z,W,iter,obj] = projectionmissingalgo2 (X,Y,lambda,numanchor,ind,lambda3,lambdaAY)
%[F,V,A,Z,iter,obj] = missingalgo_qp(X1,gt,lambda1,lambda2,ind); % X,Y,lambda,d,numanchor
% m      : the number of anchor. the size of Z is m*n.
% lambda : the hyper-parameter of regularization term.

% X      : n*di

%% initialize
maxIter = 100 ; % the number of iterations
%Weights   = opts.Weights;
m = numanchor;
numclass = length(unique(Y));
numview = length(X);
numsample = size(Y,1);
% Z = zeros(m,numsample); 
% Z(:,1:m) = eye(m);
Z = zeros(m,numsample,numview); 
AY = zeros(m,numsample);
for iiv = 1:numview
Z(:,1:m,iiv) = eye(m);
end
Y = zeros(m,numsample,numview); 
Phi = zeros(m,numsample,numview); % 张量
tol= 1e-3;
mu =1;  
max_mu = 1e10;
missingindex = constructA(ind);
for i = 1:numview
    %di = size(X{i},1); 
    A{i} = zeros(m,m); %锚点矩阵
end

for i = 1:numview
    di = size(X{i},1); 
    W{i} = zeros(di,m); %锚点矩阵
end

%alpha = ones(1,numview)/numview;
alpha = ones(1,numview);

flag = 1;
iter = 0;
%%
while flag
    iter = iter + 1;

    %% optimize Ai
    for ia = 1:numview
       %part1 = X{ia} * Z';
        part1 = W{ia}'*X{ia} * Y(:,:,ia)';
        [Unew,~,Vnew] = svd(part1,'econ');
        A{ia} = Unew*Vnew';
        A{ia}(isinf(A{ia}))=0;
        A{ia}(isnan(A{ia}))=0;
    end
  
    sX= [m,numsample,numview];
    %% optimize Z
    temp          = Y - Phi/mu;
    [vecz, ~]     = wshrinkObj(temp(:),1/mu,sX,0,3);
%   [vecz, ~]     = wshrinkObj(temp(:),Weights/mu,sX,0,2);
    Z             =  reshape(vecz, sX);









    %% optimize Y
%         C1 = zeros(:,:,numview);
%         C2 = zeros(:,:,numview);
        C1 =0;
        C2 =0;
%         for a=1:numview
%             C1=lambda*ind(:,a)'+0.5*mu * ones(1,numsample); 
%             C2=lambda* A{a}'*W{a}'*X{a}+0.5*mu*Z(:,:,a)+0.5*Phi(:,:,a);
%             
%             for ii=1:numsample
%                idx = 1:numanchor;
%                ut = C2(idx,ii)./C1(ii);
%                %Z(idx,ii) = EProjSimplex_new(ut');
%                Y(idx,ii,a)=EProjSimplex_new(ut');
%            end
%         end
    for v =1:numview
        C1 = ind(:,v)'+(mu/(2*lambda) + lambda3/lambda)*ones(1,numsample);
        C2 = (X{v})'*(W{v})*(A{v})+(lambda3/lambda)*AY'+(mu/(2*lambda))*(Z(:,:,v)+(Phi(:,:,v)./mu))';
        for ii=1:numsample
            idx = 1:numanchor;
            ut = C2(ii,idx)./C1(ii);
            Y(idx,ii,v) = EProjSimplex_new(ut');
        end
    
    end
    Y(isinf(Y))=0;
    Y(isnan(Y))=0;
    

        %% optimize AY
    %A1 = zeros(m,numsample);
    J = zeros(m,numsample); 
    for v = 1:numview
        J = J + Y(:,:,v);
    end
    J = J./numview;
%     for v = 1:numview
%         A1 = A1 + solve_l1l2(Y(:,:,v),lambdaAY/(2*lambda3));
%     end
%     AY = A1./numview;
    AY = solve_l1l2(J,lambdaAY/(2*lambda3*numview));
        
        
        
    %% optimize W
    for ib = 1:numview
        
         part2 = X{ib} * Y(:,:,ib)'*A{ib}';
         [Uneww,~,Vneww] = svd(part2,'econ');
         W{ib} = Uneww*Vneww';
    end 
    
    
    

   %%
    Phi         = Phi + mu*(Z - Y);
    chag2       = max(abs(Z(:) - Y(:)));
    
    %%
    term1 = 0;
    term2 = 0;
    term3=0;
    term4 = 0;
    for iv = 1:numview
        term1 = term1 + alpha(iv)^2 * norm(W{iv}'*X{iv} - A{iv} * (Z(:,:,iv).*repmat(missingindex{iv},m,1)),'fro')^2 ;
        term4 = term4+norm(AY-Z(:,:,iv),'fro');
        term3=term3+Z(:,:,iv);%b不用除回去
    end
    %tensor_matrix = reshape(Z, size(Z, 1), []);

%对矩阵进行奇异值分解
for t = 1:numview
[UZ, S, V] = svd(Z(:,:,t),'econ');
term2 = term2 + max(diag(S));
end
%term2 = 0;
% 计算奇异值的核范数
%nuclear_normZ = norm(S, 'fro');
column_norms = sqrt(sum(AY.^2, 1));
norm_21A = max(column_norms);
    %term3=term3/numview;
    %term2 =  norm(term3,'fro')^2;
    obj(iter) = lambda*term1+ lambda3 * term4+term2+chag2++lambdaAY*norm_21A;
     if (iter>1)&&((lambda*term1 < tol && chag2 < tol) || iter>maxIter)
          [UU,~,V]=svd(AY','econ');
        UU = UU(:,1:numclass);
         flag = 0;
     end
    
%     if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-3 || iter>maxIter || obj(iter) < 1e-10)
% %     if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-3 || iter>maxIter || obj(iter) < 1e-10)
% %        %  [UU,~,V]=svd(Z','econ');
%           [UU,~,V]=svd(term3','econ');
%         UU = UU(:,1:numclass);
%          flag = 0;
%      end
     %mu = min(mu*1.5,max_mu);
     
end
         
function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end    