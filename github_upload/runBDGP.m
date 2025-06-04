clear;
clc;

addpath(genpath('./'));
addpath(genpath('utils/'));
addpath(genpath('LRR/'));

resultdir1 = 'Results/';
if (~exist('Results', 'file'))
    mkdir('Results');
    addpath(genpath('Results/'));
end

resultdir2 = 'aResults/';
if (~exist('aResults', 'file'))
    mkdir('aResults');
    addpath(genpath('aResults/'));
end

datadir='';

dataname={'BDGP_fea'};

numdata = length(dataname); % number of the test datasets
;
numname = {'_Per0.5'};

fid=fopen('BDGP_Results_minous_best_20241123.txt','a');
fid_all=fopen('BDGP_Results_minous_all_20241123.txt','a');


for idata = 1:1 :length(dataname)
    ResBest = zeros(9, 8);
    ResStd = zeros(9, 8);
    for dataIndex = 1:1:9

        datafile = [datadir, cell2mat(dataname(idata)), cell2mat(numname(dataIndex)), '.mat'];
        load(datafile);
        %data preparation...
        gt = truelabel{1};
        cls_num = length(unique(gt));
        k= cls_num;
        tic;
       
        [X1, ind] = findindex(data, index);%index表示可观察到的
        
        time1 = toc;
        maxAcc = 0;
        %TempLambda1 = [0.1];
        TempLambda1 = [0.001  0.01 0.1  1  10];
        %TempLambda1 = [50 100 500 1000];
        TempLambda2 = [k,2*k, 3*k,5*k,7*k];
        TempLambda3 = [0.001 0.01 0.1 1 10 100];
        TempLambdaAY = (-1)*[0.001  0.01 0.1 1 10 100];
        
        ACC = zeros(length(TempLambda1),length(TempLambda2));
        NMI = zeros(length(TempLambda1), length(TempLambda2));
        Purity = zeros(length(TempLambda1), length(TempLambda2));
        %Purity = zeros(length(TempLambda1), length(TempLambda2));
        idx = 1;
        for LambdaIndex1 = 1 : length(TempLambda1)
            lambda1 = TempLambda1(LambdaIndex1);
            for LambdaIndex2 = 1 : length(TempLambda2)
                lambda2 = TempLambda2(LambdaIndex2);
                for LambdaIndex3 = 1 : length(TempLambda3)
                lambda3 = TempLambda3(LambdaIndex3);
                for LambdaIndexAY = 1 : length(TempLambdaAY)
                lambdaAY = TempLambdaAY(LambdaIndexAY);
                disp([char(dataname(idata)), char(numname(dataIndex)), '-l1=', num2str(lambda1), '-l2=', num2str(lambda2)]);
                tic;
                para.c = cls_num; % K: number of clusters
                para.k = lambda1; % m: number of nearest anchors
                [F,V,A,Z,W,iter,obj] = projectionmissingalgo2(X1,gt,lambda1,lambda2,ind,lambda3,lambdaAY); % X,Y,lambda,d,numanchor
                
                F = F ./ (eps+repmat(sqrt(sum(F .^ 2, 2)), 1, k));
                
                time2 = toc;
                stream = RandStream.getGlobalStream;
                reset(stream);
                MAXiter = 1000; % Maximum number of iterations for KMeans
                REPlic = 20; % Number of replications for KMeans
                tic;
                for rep = 1 : 20
                    pY = kmeans(F, cls_num, 'maxiter', MAXiter, 'replicates', REPlic, 'emptyaction', 'singleton');
                    res(rep, : ) = Clustering8Measure(gt, pY);
                end
                time3 = toc;
                runtime(idx) = time1 + time2 + time3/20;
                disp(['runtime:', num2str(runtime(idx))])
                idx = idx + 1;
                tempResBest(dataIndex, : ) = mean(res);
                tempResStd(dataIndex, : ) = std(res);
                ACC(LambdaIndex1, LambdaIndex2) = tempResBest(dataIndex, 1);
                NMI(LambdaIndex1, LambdaIndex2) = tempResBest(dataIndex, 2);
                Purity(LambdaIndex1, LambdaIndex2) = tempResBest(dataIndex, 3);
                Fscore(LambdaIndex1, LambdaIndex2) = tempResBest(dataIndex, 4);
                 save([resultdir1, char(dataname(idata)), char(numname(dataIndex)), '-l1=', num2str(lambda1), '-l2=', num2str(lambda2), ...
                     '-acc=', num2str(tempResBest(dataIndex,1)), '_result.mat'], 'tempResBest', 'tempResStd');
                fprintf(fid_all,'per: %f ',dataIndex);
                fprintf(fid_all,'lambda1: %f ',lambda1);
                fprintf(fid_all,'lambda2: %f ',lambda2);
                fprintf(fid_all,'lambda3: %f ',lambda3);
                fprintf(fid_all,'lambdaAY: %f ',lambdaAY);
%                fprintf(fid_all,'lambda4: %f ',lambda4);
%                fprintf(fid_all,'lambda5: %f ',lambda5);
                fprintf(fid_all,'time: %f ',runtime(idx-1));
                fprintf(fid_all,'res1: %g %g %g %g %g %g %g %g    \n',tempResBest(dataIndex, : ));
                for tempIndex = 1 : 8
                    if tempResBest(dataIndex, tempIndex) > ResBest(dataIndex, tempIndex)
                        if tempResBest(dataIndex, 1) > ResBest(dataIndex, 1)
                        l1 = lambda1;
                        l2 = lambda2;
                        l3 = lambda3;
                        lAY = lambdaAY;
                    end
                    if tempResBest(dataIndex, tempIndex) > ResBest(dataIndex, tempIndex)
                        if tempIndex == 1
                            newZ = Z;
                            newF = F;
                        end
                        ResBest(dataIndex, tempIndex) = tempResBest(dataIndex, tempIndex);
                        ResStd(dataIndex, tempIndex) = tempResStd(dataIndex, tempIndex);
                    end
                end
                end
                end
                end
            end
        end
        aRuntime = mean(runtime);
        PResBest = ResBest(dataIndex, :);
        PResStd = ResStd(dataIndex, :);
        save([resultdir2, char(dataname(idata)), char(numname(dataIndex)), 'ACC_', num2str(max(ACC(:))), '_result.mat'], 'ACC', 'NMI', 'Purity', 'aRuntime', ...
            'newZ', 'newF', 'PResBest', 'PResStd');
            fprintf(fid,'per: %f ',dataIndex);
            fprintf(fid,'lambda1: %f ',l1);
            fprintf(fid,'lambda2: %f ',l2);
            fprintf(fid,'lambda3: %f ',l3);
            fprintf(fid,'lambdaAY: %f ',lAY)
            fprintf(fid,'time: %f ',aRuntime);
            fprintf(fid,'res1: %g %g %g %g %g %g %g %g    \n',PResBest);
            fprintf(fid,'res1std: %g %g %g %g %g %g %g %g    \n',PResStd);
    end
    save();
    save([resultdir2, char(dataname(idata)), '_result.mat'], 'ResBest', 'ResStd');
end
