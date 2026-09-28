%First run 'mex transent.c' in MATLAB command window to compile base TE
%code
%Jitters source spike as in Nigam et al. 2016: "Rich-Club Organization
%[...]"
%All input units in ms.
%
% Outputs:
% 1. all_te contains TE values for all input synaptic delays in a 3D
% matrix of dimension (num_neurons,num_neurons,num_delays).
% 
% 2. ci_result contains the coincidence indexes (CI) in a
% (num_neurons,num_neurons) matrix produced by the function CIReduce.m.
%
%\rvwg2019

function [raw_te,raw_delays,jit_te,raw_ci,jit_ci] = ASDFTEr19(asdf,j_delays,i_order,j_order,windowsize,numWorkers)
    rng('shuffle')
    % Set defaults load('C:\Users\rwgarcia\Google Drive\Research\Projects\Indiana University Bloomington\c-web Similarity\old data - analyzed\02-0\asdf02-0.mat')
    if nargin < 2
        j_delays = 1;
    end

    if nargin < 3
        i_order = 1;
    end

    if nargin < 4
        j_order = 1;
    end

    if nargin < 5
        windowsize = 5;
    end

    if nargin < 6
        localCluster = parcluster('local');
        numWorkers = localCluster.NumWorkers;
    end

    num_delays = numel(j_delays);
    num_neurons = asdf{end}(1);
    
    jitstd = 5;    %standard deviation of normally-distributed jitter (in ms)
    numJit = 20;

    raw_te = zeros(num_neurons);
    raw_delays = zeros(num_neurons);
    raw_ci = zeros(num_neurons);
    
    jit_te = zeros(num_neurons,num_neurons,numJit);
    jit_ci = zeros(num_neurons,num_neurons,numJit);
    
    parpool(numWorkers)
    parfor j=1:num_neurons%source
        asdfTemp = cell(4,1);
        asdfTemp{end-1} = 1;
        asdfTemp{end}(1) = 2;
        asdfTemp{end}(2) = asdf{end}(2);

        for i=1:num_neurons%target
            if i==j
                continue
            else
                temp_raw = zeros(num_delays,1);
                
                for d=1:num_delays
                    asdfTemp([1,2]) = asdf([j,i]);  %puts source in 1 and target in 2
                    temp = transent(asdfTemp,j_delays(d),i_order,j_order);
                    temp_raw(d) = temp(2,1);
                end
                
                [raw_te(i,j),raw_delays(i,j)] = max(temp_raw);
                raw_ci(i,j) = CIReduce(temp_raw,windowsize);
                
                for s=1:numJit
                    %jitter the source neuron's spikes:
                    asdfJit = asdfTemp;
                    asdfJit{1} = unique(asdfJit{1}+...
                        round(jitstd*randn(size(asdfTemp{1}))));
                    
                    temp_jit = zeros(num_delays,1);

                    for d=1:num_delays
                        temp = transent(asdfJit,j_delays(d),i_order,j_order);
                        temp_jit(d) = temp(2,1);
                    end
                    
                    jit_te(i,j,s) = max(temp_jit);
                    jit_ci(i,j,s) = CIReduce(temp_jit,windowsize);
                end
            end
        end
        disp(j)
    end

    delete(gcp('nocreate'))
    clear localCluster ans

%     raw_ci = zeros(num_neurons);
%     if nargout>1
%         for i=1:num_neurons
%             for j=1:num_neurons
%                 if i==j
%                     continue
%                 else
%                     raw_ci(i,j) = CIReduce(raw_te(i,j,:),windowsize);
%                     
%                 end
%             end
%         end
%     end
end