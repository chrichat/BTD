function [ A,B ] = transfervariance( A,B )
%Transfer the variance from the factors in cell B to the factor A
%The norm of the columns of all the modes except from mode A are 1.
    for i=1:size(B,2)
        for ff=1:size(A,2)
            A(:,ff)=A(:,ff)*norm(B{1,i}(:,ff));
            B{1,i}(:,ff)=B{1,i}(:,ff)/norm(B{1,i}(:,ff));
        end
    end
end

