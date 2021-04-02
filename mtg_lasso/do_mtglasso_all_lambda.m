function [lambda_corrs, lambda_total, rmse_total, lambda_regs, lambda_reg_freqs, lambda_fold_regs] =  do_mtglasso_all_lambda(allY,allX,regnames,genenames, lambdas, numcvs, regweight_fname)
%% Runs multi-task group lasso to infer protein regulators for the entire module.
% 10-fold cross-validation by holding aside samples. Runs a range of lambda
% for each fold. 
% lambda_corrs: per-fold correlations, 10 per lambda
% lambda_total: concatenated correlation, 1 per lambda
% lambda_regs: avg number regulators per fold, per lambda
% lambda_reg_freqs: regulator frequencies per lambda.
% lambda_fold_regs : regualtors for each fold
%
% write regression weights for each module, regulator, fold out to
% regweight_fname

q=2;

newtargetexp=[];
for t=1:size(allY,2)
	newtargetexp=[newtargetexp;allY{t}'];
end

numsamples=size(allY{1},1);
ntasks=size(allY,2);
taskcnt=ntasks;

cvsize=floor(numsamples/numcvs);

% save regulator counts per gene across CV 
% rows are regs
% one cell per lambda.
lambda_reg_freqs=cell(size(lambdas));

% save avg number of regulators per fold, per lambda
lambda_regs=zeros(size(lambdas));

% collect module-level correlation for all lambdas
% each item in lambdas is a vector of values (size of folds)
lambda_corrs=cell(size(lambdas));    

% collect nonzero regulators for all lambdas 
% one number per fold
lambda_fold_regs=cell(size(lambdas));

% all predicted values for all test folds
lambda_predictions=cell(size(lambdas)); 

% master CC comparing all predictions to all test
lambda_total=cell(size(lambdas));
% master RMSE
rmse_total=cell(size(lambdas));

for j=1:size(lambdas,2)
    lambda_corrs{j}=zeros(numcvs,1);
    lambda_fold_regs{j}=zeros(numcvs,1);
    lambda_reg_freqs{j}=zeros(size(allX{1},2), size(allY,2));
    lambda_predictions{j}=[];
    lambda_total{j}=0;
    rmse_total{j}=0;
end

% print regression weights per fold, per lambda 
pid=fopen(regweight_fname,'w');


for cv=1:numcvs                            
	%testind=[cv];
	testbegin=((cv-1)*cvsize)+1;
	if(cv<numcvs)
		testend=cv*cvsize;
	else
		testend=numsamples;
	end
	trainind=setdiff(1:numsamples,testbegin:testend);
	%trainind=setdiff(1:numsamples,testind);
	if(isempty(trainind))
		trainind=1:numsamples;
    end
    ind=[0];
    A=[];
    Y=[];
    Atest={};
    Ytest={};
	testY=[];
	for t=1:ntasks
        Y=[Y;allY{t}(trainind,:)];
        A=[A;allX{t}(trainind,:)];
        ind=[ind size(Y,1)];
        Atest{t}=allX{t}(testbegin:testend,:);
        Ytest{t}=allY{t}(testbegin:testend,:);
		testY=[testY;allY{t}(testbegin:testend,:)'];
        %Atest{t}=allX{t}(testind,:);
        %Ytest{t}=allY{t}(testind,:);
		%testY=[testY;allY{t}(testind,:)'];
	end
    allpred=cell(size(lambdas));
    for j=1:size(allpred,2)
        allpred{j}=zeros(size(testY));
    end          

    opts=[];
    opts.init=2;
    opts.tFlag=5;
    opts.maxIter=100;   % maximum number of iterations
    opts.nFlag=0;       % without normalization
    opts.rFlag=1;       % choose max lambda and let us transit the full space
    opts.q=2;
    %opts.rho=0.0001;    % negligible impact of rho - leave empty OK?
    opts.mFlag=0;       % treating it as compositive function 
    opts.lFlag=0;       % Nemirovski's line search
    opts.ind=ind;       % group indices
    opts.fName='mtLeastR';
    tic;
    % do pathwise solution search for range of lambdas        
    X = pathSolutionLeast(A,Y,lambdas, opts);  

    % save predictions from each lambda, for each gene       
    for j=1:size(X,3)
        % LAST index is over lambdas
        x1=X(:,:,j);

        for c=1:taskcnt
            pred=Atest{c}*x1(:,c);
            allpred{j}(c,:)=pred;
            actual=Ytest{c};
            cc1=corrcoef(pred,actual);
            %fprintf('Correlations  gene %d cc1=%.2f \n',c,cc1(1,2));
        end
        
        % get module-level information        
        lampred=allpred{j};
        
        % cat 
        lambda_predictions{j}=[ lambda_predictions{j} lampred ];
        
        % flattens the two matrices into vectors
        cc=corrcoef(lampred(:),testY(:));
        %fprintf('Fold %d, Correlations modulelevel lambda=%.2f cc=%.3f\n',cv, lambdas(j), cc(1,2));

        % save value for lambda, fold
        lambda_corrs{j}(cv) =  cc(1,2);  

        % save regulators
        
        % count up number of nonzero regulators 
        
        % x1' gives coefficients for each regulator/gene pair
        % it's a matrix over regs and targets.
        for ix=1:size(x1,1)
            for jx=1:size(x1,2)
                if x1(ix,jx)~=0
                    % lambda, fold, reg, target, weight
                    %fprintf(pid, '%f\t%d\t%s\t%s\t%f\n', lambdas(j), cv, regulators.textdata{ix}, targets.textdata{jx}, x1(ix,jx));
                    %fprintf(pid, '%f\t%d\t%s\tT%d\t%f\n', lambdas(j), cv, regnames{ix}, jx, x1(ix,jx));
                    fprintf(pid, '%f\t%d\t%s\t%s\t%f\n', lambdas(j), cv, regnames{ix}, genenames{jx}, x1(ix,jx));
                end
            end
        end
        
        % which regs have *any* nonzero coeffs?
        s=sum(abs(x1'));
        num_regs=size(find(s~=0),2);          
        lambda_fold_regs{j}(cv)=num_regs;        
        
        % get coeffs that are nonzero - linear indices into x1 matrix
        nzb=find(x1~=0);           
        
        x2=zeros(size(x1));
        % set 
        x2(nzb)=1;  % mark non-zero coeffs with 1s
 
        % add to regulator counts 
        lambda_reg_freqs{j}=lambda_reg_freqs{j}+x2;       
    end   
end

% get pearson between all predictions and all test data
% get reguator frequency
for j=1:size(lambdas,2)
    
    % get pearson of all predictions and actual target data
    cc=corrcoef(lambda_predictions{j}(:), newtargetexp(:));
	lambda_total{j}=cc(1,2);
    
    % calculate RMSE per target
    rmse=sqrt(mean( (lambda_predictions{j}(:) - newtargetexp(:)).^2));
    rmse_total{j}=rmse;
        
    % get regulator frequency - *this will be returned*
    lambda_reg_freqs{j}=lambda_reg_freqs{j}/numcvs;
        
    % mean regulators per fold    
    lambda_regs(j)=mean(lambda_fold_regs{j}); %         %size(nzregs,1);
    
end
