function transmissionProbability = TE2Pij(transferEntropy, firingRate)
% transmissionProbability = transferEntropyToTransmissionProbability(transferEntropy, firingRate)
%
%    transferEntropy - (n,n) Table of transfer entropy of n neurons.
%    firingRate - (n,1) Firing rate of neurons. (firing probability in one bin)
%
% Returns:
%    transmissionProbability - (n,n) Table of transmission probability of n neurons.
%
%
% Description :
%    This is a converter from transfer entropy to transmission probability.
%    TE and TP has exact mapping.
%
%
% Example :
%    
%
% Author   : Shinya Ito
%            Indiana University
%
% Last modified on 12/18/2008

% transferEntropy = sig_TE;
% firingRate = SpontRates;
%firingRate = cellfun(@numel,asdf_raw(1:asdf_raw{end}(1)))/asdf_raw{end}(2);

% validity check
te1 = size(transferEntropy,1);
te2 = size(transferEntropy,2);
fr = size(firingRate,1);

if te1~=te2
	error('Transfer entropy must be square matrix');
elseif te1~=fr
	error('Different number of neurons between trnasfer entropy and firing rate.');
end

if nargin~=2
	error('Wrong number of arguments');
end




% number of neuron is now consistent.
n = te1;

for i=1:n
	for j=1:n
		% convert TE to TP for each neuron pair.
		% This TE from neuron j to i


		% building up TE talbe for interpolating
		fi = firingRate(i);
		fj = firingRate(j);
		c = 0:0.001:0.999; % coefficient

		b = fi - c*fj; % base firing of neuron i

		fib = 1-fi;
		fjb = 1-fj;
		cb = 1-c;
		bb = 1-b;
		bcb = 1-(b+c);

		TEtable = fjb*bb.*log2(bb./fib) + fj*bcb.*log2(bcb./fib) + fjb*b.*log2(b./fi) + fj*(b+c).*log2((b+c)./fi);
		% for debugging purpose
		% TEtable(1)=0;
		%if (rand() < 0.001)
		%	plot(c,TEtable)
		%	drawnow
		%end

		avail_ind = find((imag(TEtable)==0) .* (~(isnan(TEtable))));
		% for debugging purpose
		if length(avail_ind) < 3
			fi
			fj
			disp(avail_ind);
		else

		transmissionProbability(j,i) = interp1(TEtable(avail_ind), c(avail_ind), transferEntropy(j,i));
		end

		% avoid complex number (I'm not sure if it is good.)
		%TEtable = abs(TEtable);

		% Again this is TP from j to i
	end
end
