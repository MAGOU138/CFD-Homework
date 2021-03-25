function outputArg = f(p_,p,rho)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
global gamma;
if p_>p
    outputArg=(p_-p)/sqrt(rho*(0.5*(gamma+1)*p_+0.5*(gamma-1)*p));
else
    outputArg=-2/(gamma-1)*sqrt(gamma*p/rho)*(1-(p_/p)^(0.5-0.5/gamma));
end

end

