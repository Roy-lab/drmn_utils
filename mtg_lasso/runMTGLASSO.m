function runMTGLASSO(finname,outdir,ntasks)

% update the path to location of SLEP package
% see http://yelabs.net/software/SLEP/
addpath(genpath('SLEP_package_4.1/SLEP/'));

if ~exist(outdir)
        mkdir(outdir);
end

a=importdata(finname);
Y=a.data(:,1);
X=a.data(:,2:end);
regnames=a.textdata(1,3:end);
%regnames=a.textdata(1,2:end);
numsamples=size(Y,1)/ntasks;

genenames=a.textdata(2:numsamples+1,1);
genenames=strrep(genenames,'0.5h','');
allY=cell(1,numsamples);
allX=cell(1,numsamples);
for t=1:numsamples
	ids=t:numsamples:size(Y,1);

	allY{t} = zeroMean(Y(ids,:));
	allX{t} = zeroMean(X(ids,:));
end

seed=12345678;  % seed randomizer for reproduceability
threshold=0.05; % significance threshold for test against random baseline

lambdas=[0.1:0.1:0.9, 0.99]; % relative range for lambda
%numcvs=8;  % number of CV folds
numcvs=ntasks;  % number of CV folds
numreps=40; % random perturbations of regulator data

    
regweight_fname=sprintf('%s/regweights.tab', outdir);

resultsfname=sprintf('%s/cv_test_correlation.tab', outdir);    
plotfname=sprintf('%s/results.tab', outdir);

rmsefile=sprintf('%s/rmse.tab', outdir);

pdfname=sprintf('%s/lasso_all_lambdas.pdf', outdir);

[lambda_corrs, lambda_total, rmse_total, lambda_regs, lambda_reg_freqs, lambda_fold_regs] = do_mtglasso_all_lambda(allY,allX,regnames,genenames, lambdas, numcvs, regweight_fname);

% write and plot results 
write_per_fold_results(resultsfname, plotfname, lambdas, lambda_corrs, lambda_fold_regs);

X=cell2mat(lambda_corrs);  

fprintf('Doing randomization trials\n');

% random prefix
rand_prefix=sprintf('%s/rand', outdir);
% get random values for comparison
[rand_lambda_totals, rand_lambda_rmse, rand_lambda_reg_freqs] = eval_rand_regs(allY, allX,regnames,genenames, lambdas, numreps, numcvs, seed, @do_mtglasso_all_lambda, rand_prefix);


% lambda is success if total cc is significantly GREATER than random
% means 
%z-test: for every lambda, compare mean of real correlations to means from random runs 
successlam=[];  % indices into lambdas
labels=cell(size(lambdas));

fidrmse=fopen(rmsefile,'w');    % open RMSE file
for j=1:size(lambdas,2)
	% is corr HIGHER?
	[h,p]=ztest(lambda_total{j}, mean(rand_lambda_totals{j}), std(rand_lambda_totals{j}), 'Tail','Right');
	% is error LOWER?
	[hr,pr]=ztest(rmse_total{j}, mean(rand_lambda_rmse{j}), std(rand_lambda_rmse{j}), 'Tail','Left');
	
	fprintf('Lambda %.2f accuracy compared to random: pearson=%f, rand mu %f, rand SD %f, p=%f\n', lambdas(j), lambda_total{j}, mean(rand_lambda_totals{j}), std(rand_lambda_totals{j}), p);
	fprintf('Lambda %.2f accuracy compared to random: RMSE=%f, rand mu %f, rand SD %f, p=%f\n', lambdas(j), rmse_total{j}, mean(rand_lambda_rmse{j}), std(rand_lambda_rmse{j}), pr);
   
	fprintf(fidrmse, 'Lambda %.2f corr GREATER THAN random: pearson=%f, rand mu %f, rand SD %f, p=%f\n', lambdas(j), lambda_total{j}, mean(rand_lambda_totals{j}), std(rand_lambda_totals{j}), p);
	fprintf(fidrmse, 'Lambda %.2f err LESS THAN random: RMSE=%f, rand mu %f, rand SD %f, p=%f\n', lambdas(j), rmse_total{j}, mean(rand_lambda_rmse{j}), std(rand_lambda_rmse{j}), pr);
	 
	
	if (p < threshold || pr < threshold)
		successlam=[ successlam j ];   
		string=''; % star if corr; + if RMSE.
		if p < threshold
			string=sprintf('%s*', string);
		end
		if pr < threshold
			string=sprintf('%s+',string);
		end                           
		labels{j}=sprintf('*%.2f%s', lambda_regs(j), string);
	else
		labels{j}=sprintf('%.2f', lambda_regs(j));
	end
end   
fclose(fidrmse);    % close RMSE file

% print out regulators for successful lambdas only
if size(successlam,1) > 0
	indices=flip(successlam); % order successes biggest to smallest
	chosen=indices(1);  % start with biggest lambda
	
	% print consensus regulators for every good lambda
	for j=1:size(successlam,2)   
		fid=fopen(sprintf('%s/consensus_regs_lam%.2f.tab', outdir, lambdas(indices(j))),'w');       
		
		 % print total CC in header, with RMSE and mean regs per fold
		fprintf(fid, '# TOTAL CC\t%.2f\t%.5f\t%.5f\t%.5f\n', lambdas(indices(j)), lambda_total{indices(j)}, rmse_total{indices(j)}, mean(lambda_fold_regs{indices(j)}));
		%fprintf(fid, '# TOTAL CC\t%.2f\t%.5f\n', lambdas(indices(j)), lambda_total{indices(j)});
		fprintf(fid, '# Regulator\tFrequency\tMean_Rand_Freq\tSD_Rand_Freq\tP-val\n');
		allregs=max(lambda_reg_freqs{indices(j)}');
		nz=find(allregs>0);
		
		for s=1:size(nz,2)
			% get random mean, SD, p-value
			random_vals=rand_lambda_reg_freqs{indices(j)}(s,:);
			rand_mu=mean(random_vals);
			rand_sd=std(random_vals);
			[h, p]=ztest(allregs(nz(s)), rand_mu, rand_sd, 'Tail', 'Right');
			fprintf(fid,'%s\t%f\t%f\t%f\t%f\n',regnames{nz(s)},allregs(nz(s)), rand_mu, rand_sd, p);
		end
		fclose(fid);    
				   
	end 
end

% Plot the correlation values across lambda for real and random trials
figure;
hold on;

% plot real values
plot(cell2mat(lambda_total), '-o', 'LineWidth', 2);
%boxplot(X, 'labels',labels);
% plot random distributions; label with number of regulators
X2=cell2mat(rand_lambda_totals);
boxplot(X2, 'labels',labels);

title(sprintf('MTG-LASSO'));
xlabel('Avg. regulators per fold (* if corr. above random with p < 0.05)');
ylabel('Correlation(predicted, test)');
ylim([-1,1]);
%legend('real');
legend({'real'}, 'Position',[0.7,0.8,0.1,0.1]);
legend('boxoff');
text(4.9, 0.7, 'boxes: random');
print(gcf,'-dpdf', '-r300', sprintf('%s/mtglasso_all_lambdas.pdf', outdir));
hold off;

function zm = zeroMean(d)
m=mean(d);
m=repmat(m,size(d,1),1);
zm=d-m;
