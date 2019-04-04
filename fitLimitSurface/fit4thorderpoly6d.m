function [fittingRes] = ...
    fit4thorderpoly6d(wrench, twist, lambda, twistWeight, wrenchWeight,ifConvex)

    ifShowCode = 0;

    if norm(twist) == 0
        ifFitTwist = false;
    else
        ifFitTwist = true;
    end
%     disp('start fitting 3d/6d 4th-order polynomial');

    if size(wrench,1) == 3
        wrench = [wrench;zeros(size(wrench))];
        twist = [twist;zeros(size(twist))];
    end


    x1 = wrench(1,:);
    x2 = wrench(2,:);
    x3 = wrench(3,:);
    x4 = wrench(4,:);
    x5 = wrench(5,:);
    x6 = wrench(6,:);


    monomials4thOrder  = monomials4thorderpoly(x1,x2,x3,x4,x5,x6);

    [numDim,numSamples] = size(wrench);
    degree = 4;
    numMonomials = nchoosek(degree + numDim-1,numDim-1);
    assert(size(monomials4thOrder,2) == numMonomials );

    gradient = gradient4thorderpoly(x1,x2,x3,x4,x5,x6);

    gradientOneSample = zeros(numMonomials,numDim);
    scaling_min = eps;


%%
    cvx_begin quiet
%     cvx_begin 
        cvx_precision high
        variable Q(numDim.^2,numDim.^2) semidefinite
        variables fittedCoeff(numMonomials) wrenchFitError(numSamples) ...
            twistErrorNotNormalized(numSamples) s(numSamples) predictedTwist(numSamples,numDim)   
%         minimize(lambda * norm(fittedCoeff) + ...
%             wrenchWeight * norm(wrenchFitError) + twistWeight * norm(twistErrorNotNormalized))
        minimize(lambda * norm(fittedCoeff) + ...
            wrenchWeight * sum(wrenchFitError) + twistWeight * sum(twistErrorNotNormalized))
        
        subject to 
        for i = 1:numSamples
            norm(monomials4thOrder(i,:) * fittedCoeff - 1) <= wrenchFitError(i); % jxu: the polynomial for force
        end
        
        if ifFitTwist
            for i = 1:numSamples
            % Predicted twist, based on LS normal
                gradientOneSample(:,:) = gradient(i,:,:);
                predictedTwist(i,:) == [gradientOneSample' * fittedCoeff]'; % partial dirivative, because of velocity. Z: fitted velocity
                norm(predictedTwist(i,:) - s(i) * twist(:,i)') <= twistErrorNotNormalized(i);
%                 norm(predictedTwist(i,:) -(predictedTwist(i,:) * twist(:,i)) * twist(:,i)') <= twistErrorNotNormalized(i);
                s(i) >= scaling_min;
            end
        end
        %%
        if ifConvex
            convexConstraints4thOrder6d;
        end
    cvx_end

        % the twistError is different than the twistPreditionError, because the
    % magnitude of the predicted twist is smaller than 1
    predictedTwist = bsxfun(@rdivide, predictedTwist, sqrt(sum(predictedTwist.^2, 2)))';
    twist = bsxfun(@rdivide, twist, sqrt(sum(twist.^2)));

    twistPreditionError = sqrt(sum((predictedTwist - twist).^2));
%     [sampledWrench,shapeArray]=samplefittedlimitsurface(wrench,fittedCoeff,samplingStepSize,4);
    
    fittingRes.fittedCoeff = fittedCoeff;
    fittingRes.wrenchFitError = wrenchFitError;
    fittingRes.twistErrorNotNormalized = twistErrorNotNormalized;
    fittingRes.predictedTwist = predictedTwist;
    fittingRes.twistPreditionError = twistPreditionError;
    fittingRes.s = s;

% if ifVisualize
%       visualizedfitted6dgeometry(wrench,sampledWrench,twist);           
% %     figHandle = visualize4thorderpoly(fittedCoeff, wrench);
% end

%% Code to generate the monomials, the hessian matrix, the gradient and the convex 
if ifShowCode
    numDim = 6; % 6d 4th order polynomials
    numDegree = 4;
    syms x1 x2 x3 x4 x5 x6 z1 z2 z3 z4 z5 z6 real

    % monomials4thOrder = generatemonomials(numDegree,numVars);
    xvec=sym('x',[numDim,1],'real');
    monomials4thOrder  = loadmonomials4thorderpoly(x1,x2,x3,x4,x5,x6);

    numMonomials = size(monomials4thOrder,2);
    fittedCoeff = sym('fittedCoeff', [1 numMonomials]);  % fittedCoeff is the coeff for fitting. Final results

    % gradient of the LS => to fit the velocity and predict the velocity

    % diffX1 = diff(monomials4thOrder,x1);
    % diffX2 = diff(monomials4thOrder,x2);
    % diffX3 = diff(monomials4thOrder,x3);
    % diffX4 = diff(monomials4thOrder,x4);
    % diffX5 = diff(monomials4thOrder,x5);
    % diffX6 = diff(monomials4thOrder,x6);

    zeroArray = zeros(size(x1));
    diffX1 = [ 4.*x1.^3; 3.*x1.^2.*x2; 3.*x1.^2.*x3; 3.*x1.^2.*x4; 3.*x1.^2.*x5; 3.*x1.^2.*x6; 2.*x1.*x2.^2; 2.*x1.*x2.*x3; 2.*x1.*x2.*x4; 2.*x1.*x2.*x5; 2.*x1.*x2.*x6; 2.*x1.*x3.^2; 2.*x1.*x3.*x4; 2.*x1.*x3.*x5; 2.*x1.*x3.*x6; 2.*x1.*x4.^2; 2.*x1.*x4.*x5; 2.*x1.*x4.*x6; 2.*x1.*x5.^2; 2.*x1.*x5.*x6; 2.*x1.*x6.^2; x2.^3; x2.^2.*x3; x2.^2.*x4; x2.^2.*x5; x2.^2.*x6; x2.*x3.^2; x2.*x3.*x4; x2.*x3.*x5; x2.*x3.*x6; x2.*x4.^2; x2.*x4.*x5; x2.*x4.*x6; x2.*x5.^2; x2.*x5.*x6; x2.*x6.^2; x3.^3; x3.^2.*x4; x3.^2.*x5; x3.^2.*x6; x3.*x4.^2; x3.*x4.*x5; x3.*x4.*x6; x3.*x5.^2; x3.*x5.*x6; x3.*x6.^2; x4.^3; x4.^2.*x5; x4.^2.*x6; x4.*x5.^2; x4.*x5.*x6; x4.*x6.^2; x5.^3; x5.^2.*x6; x5.*x6.^2; x6.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX2 = [ zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; 2.*x1.^2.*x2; x1.^2.*x3; x1.^2.*x4; x1.^2.*x5; x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 3.*x1.*x2.^2; 2.*x1.*x2.*x3; 2.*x1.*x2.*x4; 2.*x1.*x2.*x5; 2.*x1.*x2.*x6; x1.*x3.^2; x1.*x3.*x4; x1.*x3.*x5; x1.*x3.*x6; x1.*x4.^2; x1.*x4.*x5; x1.*x4.*x6; x1.*x5.^2; x1.*x5.*x6; x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 4.*x2.^3; 3.*x2.^2.*x3; 3.*x2.^2.*x4; 3.*x2.^2.*x5; 3.*x2.^2.*x6; 2.*x2.*x3.^2; 2.*x2.*x3.*x4; 2.*x2.*x3.*x5; 2.*x2.*x3.*x6; 2.*x2.*x4.^2; 2.*x2.*x4.*x5; 2.*x2.*x4.*x6; 2.*x2.*x5.^2; 2.*x2.*x5.*x6; 2.*x2.*x6.^2; x3.^3; x3.^2.*x4; x3.^2.*x5; x3.^2.*x6; x3.*x4.^2; x3.*x4.*x5; x3.*x4.*x6; x3.*x5.^2; x3.*x5.*x6; x3.*x6.^2; x4.^3; x4.^2.*x5; x4.^2.*x6; x4.*x5.^2; x4.*x5.*x6; x4.*x6.^2; x5.^3; x5.^2.*x6; x5.*x6.^2; x6.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX3 = [ zeroArray; zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; x1.^2.*x2; zeroArray; zeroArray; zeroArray; 2.*x1.^2.*x3; x1.^2.*x4; x1.^2.*x5; x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x2.^2; zeroArray; zeroArray; zeroArray; 2.*x1.*x2.*x3; x1.*x2.*x4; x1.*x2.*x5; x1.*x2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 3.*x1.*x3.^2; 2.*x1.*x3.*x4; 2.*x1.*x3.*x5; 2.*x1.*x3.*x6; x1.*x4.^2; x1.*x4.*x5; x1.*x4.*x6; x1.*x5.^2; x1.*x5.*x6; x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x2.^3; zeroArray; zeroArray; zeroArray; 2.*x2.^2.*x3; x2.^2.*x4; x2.^2.*x5; x2.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 3.*x2.*x3.^2; 2.*x2.*x3.*x4; 2.*x2.*x3.*x5; 2.*x2.*x3.*x6; x2.*x4.^2; x2.*x4.*x5; x2.*x4.*x6; x2.*x5.^2; x2.*x5.*x6; x2.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 4.*x3.^3; 3.*x3.^2.*x4; 3.*x3.^2.*x5; 3.*x3.^2.*x6; 2.*x3.*x4.^2; 2.*x3.*x4.*x5; 2.*x3.*x4.*x6; 2.*x3.*x5.^2; 2.*x3.*x5.*x6; 2.*x3.*x6.^2; x4.^3; x4.^2.*x5; x4.^2.*x6; x4.*x5.^2; x4.*x5.*x6; x4.*x6.^2; x5.^3; x5.^2.*x6; x5.*x6.^2; x6.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX4 = [ zeroArray; zeroArray; zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; x1.^2.*x2; zeroArray; zeroArray; zeroArray; x1.^2.*x3; zeroArray; zeroArray; 2.*x1.^2.*x4; x1.^2.*x5; x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x2.^2; zeroArray; zeroArray; zeroArray; x1.*x2.*x3; zeroArray; zeroArray; 2.*x1.*x2.*x4; x1.*x2.*x5; x1.*x2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x3.^2; zeroArray; zeroArray; 2.*x1.*x3.*x4; x1.*x3.*x5; x1.*x3.*x6; zeroArray; zeroArray; zeroArray; 3.*x1.*x4.^2; 2.*x1.*x4.*x5; 2.*x1.*x4.*x6; x1.*x5.^2; x1.*x5.*x6; x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x2.^3; zeroArray; zeroArray; zeroArray; x2.^2.*x3; zeroArray; zeroArray; 2.*x2.^2.*x4; x2.^2.*x5; x2.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; x2.*x3.^2; zeroArray; zeroArray; 2.*x2.*x3.*x4; x2.*x3.*x5; x2.*x3.*x6; zeroArray; zeroArray; zeroArray; 3.*x2.*x4.^2; 2.*x2.*x4.*x5; 2.*x2.*x4.*x6; x2.*x5.^2; x2.*x5.*x6; x2.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x3.^3; zeroArray; zeroArray; 2.*x3.^2.*x4; x3.^2.*x5; x3.^2.*x6; zeroArray; zeroArray; zeroArray; 3.*x3.*x4.^2; 2.*x3.*x4.*x5; 2.*x3.*x4.*x6; x3.*x5.^2; x3.*x5.*x6; x3.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; 4.*x4.^3; 3.*x4.^2.*x5; 3.*x4.^2.*x6; 2.*x4.*x5.^2; 2.*x4.*x5.*x6; 2.*x4.*x6.^2; x5.^3; x5.^2.*x6; x5.*x6.^2; x6.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX5 = [ zeroArray; zeroArray; zeroArray; zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; x1.^2.*x2; zeroArray; zeroArray; zeroArray; x1.^2.*x3; zeroArray; zeroArray; x1.^2.*x4; zeroArray; 2.*x1.^2.*x5; x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x2.^2; zeroArray; zeroArray; zeroArray; x1.*x2.*x3; zeroArray; zeroArray; x1.*x2.*x4; zeroArray; 2.*x1.*x2.*x5; x1.*x2.*x6; zeroArray; zeroArray; zeroArray; x1.*x3.^2; zeroArray; zeroArray; x1.*x3.*x4; zeroArray; 2.*x1.*x3.*x5; x1.*x3.*x6; zeroArray; zeroArray; x1.*x4.^2; zeroArray; 2.*x1.*x4.*x5; x1.*x4.*x6; zeroArray; 3.*x1.*x5.^2; 2.*x1.*x5.*x6; x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; x2.^3; zeroArray; zeroArray; zeroArray; x2.^2.*x3; zeroArray; zeroArray; x2.^2.*x4; zeroArray; 2.*x2.^2.*x5; x2.^2.*x6; zeroArray; zeroArray; zeroArray; x2.*x3.^2; zeroArray; zeroArray; x2.*x3.*x4; zeroArray; 2.*x2.*x3.*x5; x2.*x3.*x6; zeroArray; zeroArray; x2.*x4.^2; zeroArray; 2.*x2.*x4.*x5; x2.*x4.*x6; zeroArray; 3.*x2.*x5.^2; 2.*x2.*x5.*x6; x2.*x6.^2; zeroArray; zeroArray; zeroArray; x3.^3; zeroArray; zeroArray; x3.^2.*x4; zeroArray; 2.*x3.^2.*x5; x3.^2.*x6; zeroArray; zeroArray; x3.*x4.^2; zeroArray; 2.*x3.*x4.*x5; x3.*x4.*x6; zeroArray; 3.*x3.*x5.^2; 2.*x3.*x5.*x6; x3.*x6.^2; zeroArray; zeroArray; x4.^3; zeroArray; 2.*x4.^2.*x5; x4.^2.*x6; zeroArray; 3.*x4.*x5.^2; 2.*x4.*x5.*x6; x4.*x6.^2; zeroArray; 4.*x5.^3; 3.*x5.^2.*x6; 2.*x5.*x6.^2; x6.^3; zeroArray]';
    diffX6 = [ zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; x1.^2.*x2; zeroArray; zeroArray; zeroArray; x1.^2.*x3; zeroArray; zeroArray; x1.^2.*x4; zeroArray; x1.^2.*x5; 2.*x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x2.^2; zeroArray; zeroArray; zeroArray; x1.*x2.*x3; zeroArray; zeroArray; x1.*x2.*x4; zeroArray; x1.*x2.*x5; 2.*x1.*x2.*x6; zeroArray; zeroArray; zeroArray; x1.*x3.^2; zeroArray; zeroArray; x1.*x3.*x4; zeroArray; x1.*x3.*x5; 2.*x1.*x3.*x6; zeroArray; zeroArray; x1.*x4.^2; zeroArray; x1.*x4.*x5; 2.*x1.*x4.*x6; zeroArray; x1.*x5.^2; 2.*x1.*x5.*x6; 3.*x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; x2.^3; zeroArray; zeroArray; zeroArray; x2.^2.*x3; zeroArray; zeroArray; x2.^2.*x4; zeroArray; x2.^2.*x5; 2.*x2.^2.*x6; zeroArray; zeroArray; zeroArray; x2.*x3.^2; zeroArray; zeroArray; x2.*x3.*x4; zeroArray; x2.*x3.*x5; 2.*x2.*x3.*x6; zeroArray; zeroArray; x2.*x4.^2; zeroArray; x2.*x4.*x5; 2.*x2.*x4.*x6; zeroArray; x2.*x5.^2; 2.*x2.*x5.*x6; 3.*x2.*x6.^2; zeroArray; zeroArray; zeroArray; x3.^3; zeroArray; zeroArray; x3.^2.*x4; zeroArray; x3.^2.*x5; 2.*x3.^2.*x6; zeroArray; zeroArray; x3.*x4.^2; zeroArray; x3.*x4.*x5; 2.*x3.*x4.*x6; zeroArray; x3.*x5.^2; 2.*x3.*x5.*x6; 3.*x3.*x6.^2; zeroArray; zeroArray; x4.^3; zeroArray; x4.^2.*x5; 2.*x4.^2.*x6; zeroArray; x4.*x5.^2; 2.*x4.*x5.*x6; 3.*x4.*x6.^2; zeroArray; x5.^3; 2.*x5.^2.*x6; 3.*x5.*x6.^2; 4.*x6.^3]';


    % this is for the convex contraints
    zvec=sym('z',[numDim,1],'real');
    % yvec = reshape(xvec*zvec',[numVars.^2,1]);
    yvec =[ x1*z1; x2*z1; x3*z1; x4*z1; x5*z1; x6*z1; x1*z2; x2*z2; x3*z2; x4*z2; x5*z2; x6*z2; x1*z3; x2*z3; x3*z3; x4*z3; x5*z3; x6*z3; x1*z4; x2*z4; x3*z4; x4*z4; x5*z4; x6*z4; x1*z5; x2*z5; x3*z5; x4*z5; x5*z5; x6*z5; x1*z6; x2*z6; x3*z6; x4*z6; x5*z6; x6*z6];


    Q = sym('Q', [numDim.^2 numDim.^2]);
    Q = tril(Q,0) + tril(Q,-1).';

    f = sum(monomials4thOrder.*fittedCoeff);
    hess = hessian(f,[x1,x2,x3,x4,x5,x6]);
   
    leftEq = simplify(expand(zvec' * hess * zvec));
    rightEq = simplify(expand(trace(yvec'*Q*yvec)));
    %% sos convex constraints
    leftEq == rightEq;
end
end

