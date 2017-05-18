function [ y ] = ifpos_damping( x1,x2 )
%ifpos_damping return x1 if x2 is positive, else zero. used for judge the
%value of the damping force
%   Detailed explanation goes here
if x2>=0&&x1>0
    y=x1;
else
    y=0;
end