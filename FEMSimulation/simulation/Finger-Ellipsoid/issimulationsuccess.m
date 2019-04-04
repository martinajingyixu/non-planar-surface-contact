function [simulationSuccess] = issimulationsuccess(simulationResultPath)

        if any(size(dir([simulationResultPath '/*.csv' ]),1))
            simulationSuccess = true;
        else
            simulationSuccess = false;
        end
end