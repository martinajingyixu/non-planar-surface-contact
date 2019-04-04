function [fittedCoeff, wrenchFitError, twistFitError, predictedTwist, s] = ...
    fit4thorderpoly3d(wrench, twist, lambda, twistWeight, wrenchWeight,ifConvex, ifVisualize, ifShowCode)

if (nargin == 7) 
    ifShowCode = 0;
end

assert(size(wrench,1) == 3,'fit 3d wrench and twist only');
assert(size(twist,1) == 3),'fit 3d wrench and twist only';
disp('start fitting 3d 4th-order polynomial');

x1 = wrench(1,:);
x2 = wrench(2,:);
x3 = wrench(3,:);

monomials4thOrder =[ x1.^4; x1.^3.*x2; x1.^3.*x3; x1.^2.*x2.^2; x1.^2.*x2.*x3; ...
    x1.^2.*x3.^2; x1.*x2.^3; x1.*x2.^2.*x3; x1.*x2.*x3.^2; x1.*x3.^3; x2.^4; x2.^3.*x3; ...
    x2.^2.*x3.^2; x2.*x3.^3; x3.^4]';

[numDim,numSamples] = size(wrench);
degree = 4;
numMonomials = nchoosek(degree + numDim-1,numDim-1);
assert(size(monomials4thOrder,2) == numMonomials );
 
zeroArray = zeros(size(x1));
diffX1 = [ 4.*x1.^3; 3.*x1.^2.*x2; 3.*x1.^2.*x3; 2.*x1.*x2.^2; 2.*x1.*x2.*x3; 2.*x1.*x3.^2; ...
    x2.^3; x2.^2.*x3; x2.*x3.^2; x3.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
diffX2 = [ zeroArray; x1.^3; zeroArray; 2.*x1.^2.*x2; x1.^2.*x3; zeroArray;...
    3.*x1.*x2.^2; 2.*x1.*x2.*x3; x1.*x3.^2; zeroArray; 4.*x2.^3; 3.*x2.^2.*x3; 2.*x2.*x3.^2; x3.^3; zeroArray]';
diffX3 = [ zeroArray; zeroArray; x1.^3; zeroArray; x1.^2.*x2; 2.*x1.^2.*x3;...
    zeroArray; x1.*x2.^2; 2.*x1.*x2.*x3; 3.*x1.*x3.^2; zeroArray; x2.^3; 2.*x2.^2.*x3; 3.*x2.*x3.^2; 4.*x3.^3]';

gradient(:,:,1) = diffX1;
gradient(:,:,2) = diffX2;
gradient(:,:,3) = diffX3;

gradientOneSample = zeros(numMonomials,numDim);
scaling_min = eps;


%%
cvx_begin quiet

    cvx_precision high
    variable Q(9,9) semidefinite
    variables fittedCoeff(numMonomials) wrenchFitError(numSamples) twistFitError(numSamples) s(numSamples) predictedTwist(numSamples,numDim)      
minimize(lambda * norm(fittedCoeff) + wrenchWeight * norm(wrenchFitError) + twistWeight * norm(twistFitError))
subject to 

    for i = 1:numSamples

        norm(monomials4thOrder(i,:) * fittedCoeff - 1) <= wrenchFitError(i); % jxu: the polynomial for force
        % Predicted twist, based on LS normal
        gradientOneSample(:,:) = gradient(i,:,:);
        predictedTwist(i,:) == [gradientOneSample' * fittedCoeff]'; % partial dirivative, because of velocity. Z: fitted velocity
        norm(predictedTwist(i,:) - s(i) * twist(:,i)') <= twistFitError(i);
        s(i) >= scaling_min;
    end
    
    if ifConvex
        convexConstraints4thOrder3d;
    end
    cvx_end

if ifVisualize
    figHandle = visualize4thorderpoly(fittedCoeff, wrench);
    visualizewrenchtwist(wrench,twist,figHandle);
end

%% Code to generate the monomials, the hessian matrix, the gradient and the convex 
if ifShowCode
    numVars = 3; % 3d 4th order polynomials
    numDegree = 4;
    syms x1 x2 x3 z1 z2 z3 yvec zvec Q fittedCoeff real

    % monomials4thOrder = generatemonomials(numDegree,numVars);
    xvec=sym('x',[numVars,1],'real');
    monomials4thOrder =[ x1^4; x1^3*x2; x1^3*x3; x1^2*x2^2; x1^2*x2*x3; ...
        x1^2*x3^2; x1*x2^3; x1*x2^2*x3; x1*x2*x3^2; x1*x3^3; x2^4; x2^3*x3; ...
        x2^2*x3^2; x2*x3^3; x3^4]';
    numMonomials = size(monomials4thOrder,2);
    fittedCoeff = sym('fittedCoeff', [1 numMonomials]);  % fittedCoeff is the coeff for fitting final results

    % gradient of the LS => to fit the velocity and predict the velocity

    % diffX1 = diff(monomials4thOrder,x1);
    % diffX2 = diff(monomials4thOrder,x2);
    % diffX3 = diff(monomials4thOrder,x3);

    zeroArray = zeros(size(x1));
    diffX1 = [ 4*x1^3; 3*x1^2*x2; 3*x1^2*x3; 2*x1*x2^2; 2*x1*x2*x3; 2*x1*x3^2; ...
        x2^3; x2^2*x3; x2*x3^2; x3^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX2 = [ zeroArray; x1^3; zeroArray; 2*x1^2*x2; x1^2*x3; zeroArray;...
        3*x1*x2^2; 2*x1*x2*x3; x1*x3^2; zeroArray; 4*x2^3; 3*x2^2*x3; 2*x2*x3^2; x3^3; zeroArray]';
    diffX3 = [ zeroArray; zeroArray; x1^3; zeroArray; x1^2*x2; 2*x1^2*x3;...
        zeroArray; x1*x2^2; 2*x1*x2*x3; 3*x1*x3^2; zeroArray; x2^3; 2*x2^2*x3; 3*x2*x3^2; 4*x3^3]';

    % this is for the convex contraints
    zvec=sym('z',[numVars,1],'real');
    % yvec = reshape(xvec*zvec',[numVars.^2,1]);
    yvec = [ x1*z1; x2*z1; x3*z1; x1*z2; x2*z2; x3*z2; x1*z3; x2*z3; x3*z3];

    Q = sym('Q', [numVars.^2 numVars.^2]);
    Q = tril(Q,0) + tril(Q,-1).';

    f = sum(monomials4thOrder.*fittedCoeff);
    hess = hessian(f,[x1,x2,x3]);
   
    leftEq = simplify(expand(zvec' * hess * zvec));
    rightEq = simplify(expand(trace(yvec'*Q*yvec)));
    %% sos convex constraints
    leftEq == rightEq;
end


end

