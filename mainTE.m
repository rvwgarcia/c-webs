% Compile first: mex transent1dir.c and transent.c
% \rvwg2019

% addpath ./TEpackage_old;
% addpath ./TEpackage

% Run me from '...\Research\Projects\Indiana University Bloomington\SOC by c-webs'


%% 1. Calculate TE for each delay:
options = 'load';
delays = 1:16;
dataset = '18_06-0';
dir = '..\c-web Similarity\old data - analyzed\';
load(strcat(dir,dataset,'\asdf',dataset(4:end)),'asdf_raw')

asdf = asdf_raw;
N = asdf{end}(1);

if strcmp(options,'zip')
    asdf_raw(1:N) = cellfun(@(x) cast(x,'int32'),asdf_raw(1:N),'un',0); %zip the data
    save(strcat(dir,dataset,'\asdf',dataset(4:end)),'asdf_raw')
    clear asdf_raw
elseif strcmp(options,'unzip')
    asdf(1:N) = cellfun(@(x) cast(x,'double'),asdf(1:N),'un',0);    %unzip the data
elseif strcmp(options,'load')
    asdf = asdf_raw;
    clear asdf_raw
end

%% Run TE:
tic
[raw_te,raw_delays,jit_te,raw_ci,jit_ci] = ASDFTEr19(asdf,delays);
runtime = toc;

fname = ['tempj20r19setb',dataset];
save(fname,'raw_te','raw_delays','jit_te','raw_ci','jit_ci','runtime') %jit_te may be too large to save (>2gb)

%% Determine significant connections. (1) First, plot CI against log(TEmax)
% for all neuron pairs from actual data.

figure
imagesc(raw_te) %only zeros are the N=180 along the diagonal
N = size(raw_te,1);

figure
scatter(log10(reshape(raw_te,[N^2,1])),reshape(raw_ci,[N^2,1]),'.')

% (2) Plot CI vs log(TEmax) for jittered data.
numJit = size(jit_te,3);
for jitter=1:numJit
    hold on
    scatter(log10(reshape(jit_te(:,:,jitter),[N^2,1])),reshape(jit_ci(:,:,jitter),[N^2,1]),'ko')
end

set(gca,'FontSize',12)
title(strcat('Data set:',[' ' dataset(1:2) '\' dataset(3:end)]),...
    'interpreter','latex','FontSize',12)
xlabel('$\log_{10}($TE$^*)$','interpreter','latex','FontSize',12)
ylabel('Coincidence index, CI','interpreter','latex','FontSize',12)
legend('Raw data','Jittered')
savefig(strcat('sigTE-CI',dataset))

% (3) Combine the raw and jittered plots to construct a filter to indicate
% the ratio of jittered to all dots in each of the 25x25 bins.

maxCI = max([max(jit_ci(jit_ci>0)) max(raw_ci(raw_ci>0))]);
minCI = min([min(jit_ci(jit_ci>0)) min(raw_ci(raw_ci>0))]);

maxTE = max([max(jit_te(jit_te>0)) max(raw_te(raw_te>0))]);
minTE = min([min(jit_te(jit_te>0)) min(raw_te(raw_te>0))]);

if numel(jit_ci(jit_ci>0))~=numel(jit_te(jit_te>0))
    warning('The number of TE$>0$ values is not equal to the number of CI$>0$ values.')
    n = numel(jit_te(jit_te>0));
    m = numel(jit_ci(jit_ci>0));
    disp(n)
    disp(m)
else
    n = numel(jit_te(jit_te>0));
end

numBins = ceil((2*n)^(1/3));  %cf. Terrell & Scott 1985
numBins = 200;

CI = linspace(minCI,maxCI,numBins+1);
TE = linspace(minTE,maxTE,numBins+1);
E = zeros(numBins);

tic
for t=1:numBins
    for c=1:numBins
        tempJit = intersect(find(jit_te>=TE(t) & jit_te<TE(t+1)),...
            find(jit_ci>=CI(c) & jit_ci<CI(c+1)));
        tempRaw = intersect(find(raw_te>=TE(t) & raw_te<TE(t+1)),...
            find(raw_ci>=CI(c) & raw_ci<CI(c+1)));
        
        Njit = numel(tempJit);
        Nraw = numel(tempRaw);
        E(t,c) = Njit/(Nraw+Njit);
    end
end
toc

E(isnan(E)) = 1;

figure
imagesc(log10(TE),CI,log10(E'))
set(gca,'Ydir','normal')
set(gca,'FontSize',12)
title(strcat('Data set:',[' ' dataset(1:2) '\' dataset(3:end)]),...
    'interpreter','latex','FontSize',12)
xlabel('$\log_{10}($TE$^*)$','interpreter','latex','FontSize',12)
ylabel('Coincidence index, CI','interpreter','latex','FontSize',12)
colormap('gray');
h = colorbar;
ylabel(h,'$\log_{10}($error rate$)$','interpreter','latex','FontSize',12)
savefig(strcat('errorBinsTE-CI',dataset))

% (4) Real connections are presumed when the error rate is below 0.03
% (Shimono & Beggs 2015)

Eth = 0.03;
[TEind,CIind] = find(E<Eth);

sig_te = zeros(N);

for ind=1:numel(TEind)
    t = TEind(ind);
    c = CIind(ind);
    
    tempRaw = intersect(find(raw_te>=TE(t) & raw_te<TE(t+1)),...
        find(raw_ci>=CI(c) & raw_ci<CI(c+1)));
    
    for k=1:numel(tempRaw)
        [i,j] = ind2sub([N,N],tempRaw(k));
        sig_te(i,j) = raw_te(i,j);
    end
end

figure
imagesc(sig_te)

%%
%%\section{1} Using the data just produced (raw_te,jit_te,ci_result), this
%first algorithm eliminates self-connections prior to jittering:
% raw_te = permute(raw_te,[2,1,3]); %files in 'temp1E2r19' need this, but
% this error is now fixed in the latest version of ASDFTEr19
% jit_te = permute(jit_te,[2,1,3,4]);
% ci_result = ci_result';
[max_te,del_te] = max(raw_te,[],3); %finds the largest TE value along the delays

%set target neuron
a = 1;
N = size(max_te,1);

figure
plot(max_te(a,:))   %plot against other input TE values
hold on
plot(max_te(a,a)*ones(N,1),'k--')
xlabel('Source neuron, $j$','interpreter','latex')
ylabel('TE value, $T_{aj}$','interpreter','latex')
legend(strcat('TE for target neuron $a=',num2str(a),'$'),'Self-connection amplitude, $T_{aa}$')
xlim([1,N])

%subtract-off the self-connection amplitude
Ath = zeros(size(max_te));

for i=1:N
    Ath(i,:) = max_te(i,:)-max_te(i,i);
end

Ath(Ath<0) = 0; %but this eliminates some of the previously-identified
                %connections in sig_te:
Bcorr = Ath.*sig_te;
CorrAmp = numel(find(Bcorr))/N^2*100; %There is only ~1% overlap...

%%\section{2} This second algorithm utilizes the jit_te prior to
%elimination of self-connections:
raw_te_sig = zeros(N,N,numel(delays));
tic
for i=1:N
    for j=1:N
        for d=1:numel(delays)
            temp = squeeze(jit_te(i,j,d,:));
            raw_te_sig(i,j,d) = raw_te(i,j,d)-mean(temp);
            
            if raw_te_sig(i,j,d)<std(temp)
                raw_te_sig(i,j,d) = 0;
            else
                continue
            end
        end
    end
    disp(i)
end
toc

[max_te2,del_te2] = max(raw_te_sig,[],3);

figure
plot(max_te2(a,:))   %plot against other input TE values
hold on
plot(max_te2(a,a)*ones(N,1),'k--')
xlabel('Source neuron, $j$','interpreter','latex')
ylabel('TE value, $T_{aj}$','interpreter','latex')
legend(strcat('TE for target neuron $a=',num2str(a),'$'),'Self-connection amplitude, $T_{aa}$')
xlim([1,N])

%Now subtract-off the self-connection:
Ath2 = zeros(size(max_te2));

for i=1:N
    Ath2(i,:) = max_te2(i,:)-max_te2(i,i);
end

Ath2(Ath2<0) = 0; %but this eliminates some of the previously-identified
                %connections in sig_te:
Bcorr2 = Ath2.*sig_te;
CorrAmp2 = numel(find(Bcorr2))/N^2*100; %There is only ~3% overlap...

%connectivity fractions:
disp(numel(find(Ath))/N^2*100) %~18.0%
disp(numel(find(Ath2))/N^2*100) %~64.5%
disp(numel(find(sig_te))/N^2*100) %~3.7%, which falls within the range of
                                  %expected connectivity density of an
                                  %effective network of Izhikevich neurons
                                  %(Nigam et al. 2016)

save('sigTEresults_temp1E2r19','Ath','Ath2','sig_te')

ind = find(Ath2); disp(numel(ind)/N^2*100)
TEmax = Ath2(ind);
CI = ci_result(ind);

figure
scatter(TEmax,CI)
xlabel('TE value, $T_{ij}$','interpreter','latex')
ylabel('CI value, $C_{ij}$','interpreter','latex')

ind = find(Ath); disp(numel(ind)/N^2*100)
TEmax = Ath(ind);
CI = ci_result(ind);
hold on
scatter(TEmax,CI,'s')

ind = find(sig_te); disp(numel(ind)/N^2*100)
TEmax = sig_te(ind);
CI = ci_result(ind);
hold on
scatter(TEmax,CI,'x')
xlim([0,2E-5])
%Now we must develop a way to determine the boundary between the true
%connections and the false positives

% str = ['mkdir ./data/',num2str(data_name),'_graph/;'];   eval(str);
% str = ['save ./data/',num2str(data_name),'_graph/TEdelays TEdelays;'];   eval(str);
%clear TEdelays_rand_all TEdelays0_rand

%% Shuffle the spike times to test for TE value signiifice
% clear TEdelays_rand_all TEstds_rand_all ;

numShuf = 1E3;
time_ax = numel(delays);
TEdelays_rand_all = zeros(N,N,time_ax);

tic
for k=1:numShuf
    asdf_rand = asdf;
    for kk=1:N
        for jj=1:size(asdf{kk},2)-5
            jitter_bin = randperm(19,1)-10;
            if asdf_rand{kk}(jj)>=2
                asdf_rand{kk}(jj) = asdf{kk}(jj) + jitter_bin;
                asdf_rand{kk}(jj+1) = asdf{kk}(jj+1) - jitter_bin;
            end
        end
        asdf_rand{kk} = sort(asdf_rand{kk});
        
        asdf_mod = cell(size(asdf));
        asdf_mod{1} = asdf_rand{kk};
        
        ind = setdiff(1:N,kk);
        temp = ASDFSubsample(asdf,ind);
        temp{end}(1) = temp{end}(1)+1;
        asdf_mod(2:end) = temp;
        
        TEdelays_rand = ASDFTE_parallel_mod(asdf_mod,delays);
%         TEdelays_rand = ASDFTE(asdf_mod,delays);
%         clear asdf_mod
        TEdelays_rand_all(kk,:,:) = TEdelays_rand_all(kk,:,:)+TEdelays_rand;
%         clear TEdelays_rand
    end
    disp(k)
%     str = ['save ./data/',num2str(data_name),'_graph/TEdelays_rand_all',num2str(k),' TEdelays_rand_all;']; eval(str);
end
% clear TEdelays_rand_all
% save('temp','TEdelays_rand_all')
toc
runtime = toc;

TEdelays_rand_all = TEdelays_rand_all/numShuf;

[sig_te,sig_del] = max(all_te-TEdelays_rand_all,[],3);
sig_te(sig_te<0) = 0;
sig_del(sig_te==0) = 0;

delete(gcp('nocreate'))

save('temp1E3')
