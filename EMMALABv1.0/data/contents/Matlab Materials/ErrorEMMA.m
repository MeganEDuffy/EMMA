classdef ErrorEMMA < matlab.mixin.SetGet & matlab.mixin.Copyable
    %Class to calculate jackknifed endmember values.
    %   This class is only used if the user wishes to perform an error
    %   analysis on their endmember datasets. This is done by jackknifing
    %   each endmember groups by their mean, thus calculating the error
    %   within the endmember datasets themselves. This can generally be
    %   used to suggest whether or not the group of endmembers are
    %   representative of each other. 

    properties
        Data DataEMMA
        Analyze AnalyzeEMMA
        EMDataNum (1, :) cell %Each cell contains the numeric data for each selected EM. 
        EMZData (1, :) cell %Each cell conatins the z-scored numeric data for the jackknifed EMs.
        JackknifeEM (1, :) cell %
        Combine double %Matrix of all possible combinations of endmembers
        MeanPercents double
        StdPercents double
        Time double %The amount of time it took for the whole sequence to run
        MeanData (1, :) cell %Mean endmember contributions separated by location
        StdData (1, :) cell %Std dev endmember error separated by location 
    end

    methods
        function obj = ErrorEMMA(deData, aData)
            %ErrorEMMA accepts both a DataEMMA and AnalyzeEMMA value
            arguments
                deData DataEMMA
                aData AnalyzeEMMA
            end

            obj.Data = deData;
            obj.Analyze = aData;
            obj = Jackknife(obj);
            obj = RepMix(obj);
        end

        function obj = Jackknife(obj)
            %Jackknife each group of the selected EM and store in cell
            %array.
            obj.Data.FullNumericLoc = cell(1, length(obj.Data.SelectedEMNum));
            
            for i = 1:length(obj.Data.SelectedEMNum)
                obj.JackknifeEM{i} = jackknife(@mean, obj.Data.SelectedEMNum{i});
            end
            
            %Combine all the data back together to z-score.
            mNumAllEM = obj.JackknifeEM{1};

            for i = 2:length(obj.Data.SelectedEMNum)
                mNumAllEM = [mNumAllEM; obj.JackknifeEM{i}];
            end
            
            %Z-score the jackknifed endmember values using the mean and
            %standard deviation of the mixture data.
            mEndmemberJackZ = zeros(size(mNumAllEM));

            for i = 1:length(obj.Data.mu)
                mEndmemberJackZ(:, i) = (mNumAllEM(:, i) - obj.Data.mu(i))/obj.Data.sigma(i);
            end

            %Separate data back out again.
            %Store numerical endmember data into cell by endmember type
            iRowOn = 0;
            cEMZData = cell(1, length(obj.Data.SelectedEMNum));

            for i = 1:length(obj.Data.SelectedEMNum)
                iRows = size(obj.Data.SelectedEMNum{i}, 1);
                cEMZData{i} = mEndmemberJackZ(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end
            
            obj.EMZData = cEMZData;

            %Separate out the sizes of the jackknifed endmember values.
            vSize = cell(1, length(obj.Data.SelectedEMNum));
            vSize{1} = 1:size(obj.Data.SelectedEMNum{1}, 1);
            
            for i = 2:length(obj.Data.SelectedEMNum)
                vSize{i} = 1:size(obj.Data.SelectedEMNum{i}, 1); %This might give trouble in debugging lol. 
            end

            %Create all possible combinations of the selected jackknifed
            %endmembers
            vCombine = vSize{1};

            for i = 2:length(obj.Data.SelectedEMNum)
                vCombine = combvec(vCombine, vSize{i});
            end

            obj.Combine = vCombine;
        end

        function obj = RepMix(obj)
            %Time how long it takes for the process
            tic

            %Calculate the percentages for the different EM combinations.
            mA = obj.Data.AvgEMSelectedZ';
            mB = obj.Data.FullMixSelectedZ'; 
            vEMSel = length(obj.Data.SelectedEMNum);
            vEMData = obj.EMZData;
            vCombine = obj.Combine;

            iEndMbr = size(mA, 2); 
            iSamples = size(mB, 2);

            mPercents = zeros(iSamples, iEndMbr);

            mMeanPercents = zeros(iEndMbr, iSamples);
            mStdPercents = zeros(iEndMbr, iSamples);

            for i = 1:iSamples
                vB = mB(:, i);
                mPercents(i, :) = obj.Analyze.MixCalc(mA, vB);
            
                iNumPerm = size(vCombine, 2); 
            
                % Create a zeros matrix for the calculated values
                mPercentsIter = zeros(size(vCombine')); 

                parfor n = 1:iNumPerm % this is really big, so consider doing parfor here instead of for

                    vJackIter = vEMData{1}(vCombine(1, n), :);

                    for y = 2:vEMSel
                    vJackNext = vEMData{y}(vCombine(y, n), :);
                    vJackIter = [vJackIter; vJackNext];
                    end
            
                    % Run the mixing calculation and store the results.
                    mPercentsIter(n, :) = obj.Analyze.MixCalc(vJackIter', vB);
                    % Put the percents calculated in column i of mPercentsInter
                
                end
            
                mMeanPercents(:, i) = mean(mPercentsIter', 2);
                mStdPercents(:, i) = std(mPercentsIter', 0, 2);
            
            end 
            
            mMeanPercents = mMeanPercents';
            mStdPercents = mStdPercents';

            obj.MeanPercents = mMeanPercents;
            obj.StdPercents = mStdPercents;

            %End of time elapsed
            toc
            obj.Time = toc;


            %Separate out the standard deviation and mean percents by
            %sample location.
            iRowOn = 0;
            obj.MeanData = cell(1, length(obj.Data.SelectedLocData));

            for i = 1:length(obj.Data.SelectedLocData)
                iRows = size(obj.Data.SelectedLocData{i}, 1);
                obj.MeanData{i} = obj.MeanPercents(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end

            iRowOn = 0;
            obj.StdData = cell(1, length(obj.Data.SelectedLocData));

            for i = 1:length(obj.Data.SelectedLocData)
                iRows = size(obj.Data.SelectedLocData{i}, 1);
                obj.StdData{i} = obj.StdPercents(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end
            
        end
    end

    methods (Static)
        
    end
end