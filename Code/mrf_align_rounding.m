function [X_int] = mrf_align_rounding(Data, X, y, lambda)

if 1
    dimX = size(Data.W_X, 1)*size(Data.W_X, 2);
    numCons = length(Data.idsF);
    Data.F = sparse(1:numCons,...
        Data.idsF,...
        ones(1, numCons),...
        numCons, dimX);
    Data.Ft_full = Data.F';
    % edge based rounding
    [numS, numT] = size(X);
    scoreY = Data.Ft_full*Data.E*y;
    score = lambda*scoreY + reshape(X, [dimX, 1]);
    [s, order] = sort(-score');
    
    flags_S = zeros(1, numS);
    flags_T = zeros(1, numT);
    
    X_int = zeros(numS, numT);
    count = 0;
    for i = 1:length(order)
        if count == min(numS, numT)
            continue;
        end
        corId = order(i);
        tId = floor((corId-1)/numS) + 1;
        sId = corId - (tId-1)*numS;
        if flags_S(sId) == 0 && flags_T(tId) == 0
            X_int(sId, tId) = 1;
            flags_S(sId) = 1;
            flags_T(tId) = 1;
            count = count + 1;
        end
    end
else
    
end



