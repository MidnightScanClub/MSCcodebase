function [z] = FisherTransform(r)
%This function performs an r to z (Fisher) transformation

r(r > .999) = .999;%logical indexing: the test makes an array of size size(r), type logical, with 1 where r == 1, 0 otherwise
%when using a logical array of matching dimension as index, it performs a
%specialized indexing, returning a vector that refers to only those
%elements of the original matrix

%{
if r == 1
    r = .999;
end
%}

z = 0.5 * log((1+r)./(1-r));

end

