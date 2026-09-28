% Returns weighted adjacency matrix at maximum dynamic range using the
% Perron-Frobenius eigenvalue, as per Larremore et al. Phys Rev Lett 106,
% 058101 (2011). The output 'red' indicates a reducible matrix (1) or
% irreducible matrix (0).
% 
% Rashid Williams-Garcia 2015

function [A,PFeig0,red] = PerronFrobEigvl(A)
    PFeig0 = max(abs(eig(A))); %original Perron-Frobenius eigenvalue
	A = A/PFeig0;
    
    N = length(A);
    B = (eye(N)+A)^(N-1);
    Btest = B>0;
    if numel(Btest(Btest==1))==N^2
        red = 0;
        display('The adjacency matrix is irreducible.')
    else
        red = 1;
        display('The adjacency matrix is reducible.')
    end
end