classdef AnalyzeEMMA < matlab.mixin.SetGet & matlab.mixin.Copyable
    %Calculations used to generate endmember fractional contributions
    %   AnalyzeEMMA reads in a DataEMMA value, which is then used to
    %   calculate the fractional contributions of the users selected
    %   endmembers. The MATLAB function 'fmincon' is used to generate these
    %   contributions. Additionally, predicted values based on the
    %   generated fractional contributions are calculated to assess how
    %   well the model predicted endmember contributions. 

    properties
        Data DataEMMA
        EMData double
        MixData double
        MixturePercentages (:, :) double %Matrix of percent contributions for each endmember over time
        FracData (1, :) cell %Cell array where each cell is the fractional contributions by location
        PredictedVal double %Predicted values for all numeric data
        PredictionLoc (1, :) cell %Cell array of prediction values to allow for easy plotting in app.
    end

    methods
        function obj = AnalyzeEMMA(deData)
            %AnalyzeEMMA accepts a DataEMMA value. 
             
            arguments
                deData DataEMMA
            end

            obj.Data = deData;

            obj.EMData = obj.Data.AvgEMSelectedZ;
            obj.MixData = obj.Data.FullMixSelectedZ; 

            obj = Optimization(obj);
            obj = PredictCalc(obj);
        end

        function obj = set(obj, sIndex, value)
            switch sIndex
                case 'EMData'
                    obj.EMData = value;

                case 'MixData'
                    obj.MixData = value;
                    obj = Optimization(obj);
                    obj = PredictCalc(obj);
            end
        end

        function obj = Optimization(obj)
            iEndMbr = size(obj.EMData, 1); 
            iSamples = size(obj.MixData, 1);
            
            obj.MixturePercentages = zeros(iSamples, iEndMbr);
            mB = obj.MixData';
            mA = obj.EMData';
            %If alternative option to solve through direct solution, it
            %would go here with an 'if' statement

                for i = 1:iSamples
                    vB = mB(:, i);
                    obj.MixturePercentages(i, :) = MixCalc(obj, mA, vB);
                end
            %end
            %Store fractional contributions of each endmember into a cell
            %array for easy data organization. 
            iRowOn = 0;
            obj.FracData = cell(1, length(obj.Data.SelectedLocData));

            for i = 1:length(obj.Data.SelectedLocData)
                iRows = size(obj.Data.SelectedLocData{i}, 1);
                obj.FracData{i} = obj.MixturePercentages(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end
        end 

        function vMix = MixCalc(obj, mA, vB)
            %Initial guesses for the percent contributed from each endmember. x0 must
            %have the same number of rows as mA has columns.
            % x0 = [0.2 0.2 0.2 0.2 0.2]';
            %How many rows in x0? 
            iRows = size(mA, 2);
            numFill = 1/iRows; 
            x0 = numFill * ones(iRows, 1);
            
            %Call fmincon to minimize.
            % Start loop here for creating array of percentages.
            Aeq = ones(1, iRows);
            
            [vMix, ~] = fmincon(@(x) obj.ObFunc(x, mA, vB), x0, [], [], Aeq, 1,... % change ones vector based on # of endmembers
                zeros(size(x0)), ones(size(x0)));
        end

        function obj = PredictCalc(obj)
            %Calculate the predicted values for a selected tracer.
            %Convert percentages back to z-scored values
            mZEMData = obj.Data.AvgEMSelectedZ;
            iTracers = size(mZEMData, 2);
            iSamples = size(obj.MixData, 1);
            mCalculatedZ = zeros(iSamples, iTracers);
            mFractions = obj.MixturePercentages;
            
                for i = 1:iSamples
                    mCalculatedZ(i, :) = mFractions(i, :) * mZEMData;
                end
            
            %Unzscore the predicted values. This is done by multipling the original
            %values by the standard deviation (sigma) and adding the mean (mu). 
            mAdjusted = zeros(size(mCalculatedZ));
            
                for i = 1:iTracers
                    mAdjusted(:, i) = mCalculatedZ(:, i) * obj.Data.sigma(i);
                    mAdjusted(:, i) = mAdjusted(:, i) + obj.Data.mu(i);
                end
            
            obj.PredictedVal = mAdjusted;

            %Store numerical prediction data into cell by location.
            iRowOn = 0;
            cPrData = cell(1, length(obj.Data.MixturePCA));

            for i = 1:length(obj.Data.MixturePCA)
                iRows = size(obj.Data.MixturePCA{i}, 1);
                cPrData{i} = mAdjusted(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end
 
            obj.PredictionLoc = cPrData;
        end
    end

    methods (Static)
        function dSqErr = ObFunc(vX, mA, vB)
            %Sum of squared errors.  
            vBcalc = mA * vX;
            dSqErr = sum((vBcalc - vB).^2); 
        end
    end
end