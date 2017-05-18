function [ b ] = ifpos( a )
%IFPOS tests if the input arg is positive. If positive, return the same
%value, else return zero.
b = (a+abs(a))/2;


end

