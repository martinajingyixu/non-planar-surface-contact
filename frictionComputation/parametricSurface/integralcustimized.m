function [integratedValue] = integralcustimized(integrand,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)
% check if the integrand is about 0,1,or 2 variables.
% the name of the variable can only be u or v.

syms u v
variables = symvar(integrand);

if numel(variables) == 2

    integratedValue = integral2(matlabFunction(integrand)...
                ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);
elseif numel(variables) == 0
%         integratedValue = 0;
    uint = integral(@(u)(double(integrand)+0*u),lowerBoundU,upperBoundU);
    integratedValue = integral(@(v)(uint)+0*v,lowerBoundV,upperBoundV);
else
    if strcmp(char(variables(1)),'u')
        uint = integral(matlabFunction(integrand)...
            ,lowerBoundU,upperBoundU);
        integratedValue = integral(@(v)(uint+0*v),lowerBoundV,upperBoundV);
    elseif strcmp(char(variables(1)),'v')
        vint = integral(matlabFunction(integrand)...
            ,lowerBoundV,upperBoundV);
        integratedValue = integral(@(u)(vint+0*u),lowerBoundU,upperBoundU);

    else
        assert("ERROR! Can only have variables u or v");
    end
end
    
end

