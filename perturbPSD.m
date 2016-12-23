function [ Sigma ] = perturbPSD( Sigma, minEig )
% With inputs of a symmetric matrix Sigma and its minimum eigen value,
% perturn Sigma by a bit to return a symmetric positive-semidefinite matrix.
% Sources:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/103174 (with
% slight change).
r = size(Sigma,1);
Sigma = Sigma - (Sigma<0).*Sigma; % Make sure all entries are positive
if (minEig < 0)
    absMinEig = abs(minEig);
    triag = triu(Sigma);
    triag = triag + (absMinEig+0.000001)*eye(r);
    triag = triag./(1+absMinEig+0.000001);
    % copy the upper triagle to lower to make sure symmetric
    Sigma = triag + triu(triag,1)'; 
end
end

