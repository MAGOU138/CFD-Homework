function [outputArg,c] = solve_F(du,pL,pR,rhoL,rhoR)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

p0=1.e-8;
dp=1.e-8;

delta=1;
count=0;
while delta>1.e-8
    count=count+1;
    k=((f(p0+dp,pL,rhoL)+f(p0+dp,pR,rhoR))-(f(p0,pL,rhoL)+f(p0,pR,rhoR)))/dp;
    p1=p0-(f(p0,pL,rhoL)+f(p0,pR,rhoR)-du)/k;
    delta=abs(p1-p0);
    p0=p1;
end

outputArg=p0;
c=count;

end

