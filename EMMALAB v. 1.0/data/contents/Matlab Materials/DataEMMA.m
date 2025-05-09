classdef DataEMMA < matlab.mixin.SetGet & matlab.mixin.Copyable
    %DataEMMA generates all necessary variables to start the mixing
    %calculation. 
    %   This class generates properties that are required for the mixing
    %   analysis. In order to allow for ease with data organization and
    %   user interaction within the app EMMA, a range of properties are
    %   generated. The standardization of both mixture and endmembers is
    %   calculated here, along with the PCA on the mixture and endmember
    %   data.

    properties
        FullPath char %Path to data file
        DataAll table %Table of all data, cleaned
        MixTable table %Table of mixture data, separated from DataAll
        EMTable table %Table of endmember data, separated from DataAll
        Locations (1, :) cell %List of different mixture locations 
        NumLoc int16 %Number of locations in the mixture data
        TypeEM (1, :) cell %List of different endmember types
        NumEM int16 %Number of endmembers 
        MixtureLoc (1, :) cell %Each cell has a table with the data for each location in the mixture data
        Endmembers (1, :) cell %Each cell has a table with the data for each endmember
        Tracers (1, :) cell %Each cell is the tracer defined by the column names for the elements or isotopes in DataAll
        StartDate datetime %Selected start date to view data (default is the earliest date in mixture dataset)
        EndDate datetime %Selected end date to view data (default is the earliest date in mixture dataset)
        SelectedEM (1, :) logical %Each cell is a logical vector generated from the users selection on the table checkbox. 
        SelectedLoc (1, :) logical %Each cell is the name of a single location selected
        SelectedLocData (1, :) cell %Each cell is the selected location with the specified date range.
        SelectedMixDataAll table %Full table of selected data by location. 
        SelectedEMData (1, :) cell %Cells containing the endmember data for the selected endmembers in the datasets.
        FullMixNumeric double %Contains the full numeric data for the SELECTED locations in a single matrix. 
        FullNumericLoc (1, :) cell %Contains the full numeric data for the selected locations, where each cell is a different location.
        SelectedEMNum (1, :) cell %Contains the full numeric data for the selected endmembers, where each cell is a different endmember.
        CorrMatrix double %Matrix of values
        SelectedTracers (1, :) logical %Each cell is a logical vector generated from the users selection on the table checkbox. 
        MixturePCA (1, :) cell %Should this be separated with each cell being a different location as well?
        VariablesPCA table %table of PCA scores with PC titles
        MixPCAData struct %Cell array containing, in order, the coefficient matrix, the percent explained, and the ____ for the mixture data. 
        EMPCA (1, :) cell %Same questions as above.
        EMScoreMatrix table %Endmember data but as a matrix of all values
        AvgEMScore table %Averaged endmember dataset in matrix form.
        AvgEMPCA (1, :) cell %Each cell is an averaged endmember with each column being a selected tracer
        FullMixSelectedZ double %Save properties for full selected dataset and averaged selected
        AvgEMSelectedZ double %endmembers as zscored values
        mu double
        sigma double
        NumSelectLoc (1, :) cell %what's the difference between this and SelectedLocData?
    end

    methods
        function obj = DataEMMA(sFullPath)
            %DataEMMA accepts a table of data
            %   The input into DataEMMA is the table of data, arranged
            %   according to the instructions for the Microsoft Excel
            %   template.

            arguments 
                sFullPath char = ''
            end

            %If a path was passed, load that file. If not, open a dialog
            %box to find the file to load. 
            if isempty(sFullPath) 
                [sFile, sPath] = uigetfile('*.xlsx', 'Select Data File');
                if sFile == 0
                    return;
                end

                sFullPath = fullfile(sPath, sFile); 
                
            end

            obj.FullPath = sFullPath;

            %Load the data. 
            obj = LoadData(obj); 
            obj = FindCat(obj);
            obj = SeparateData(obj);
            obj = CalcPCA(obj);

        end

        function obj = set(obj, sIndex, value)
            switch sIndex
                case 'StartDate'
                    obj.StartDate = value;
                    obj = SeparateData(obj);
                    obj = CalcPCA(obj);

                case 'EndDate'
                    obj.EndDate = value;
                    obj = SeparateData(obj);
                    obj = CalcPCA(obj);

                case 'SelectedLoc'
                    obj.SelectedLoc = value;
                    obj = SeparateData(obj);
                    obj = CalcPCA(obj);

                case 'SelectedEM'
                    obj.SelectedEM = value;
                    obj = CalcPCA(obj);

                case 'SelectedTracers'
                    obj.SelectedTracers = value;
                    obj = SeparateData(obj);
                    obj = CalcPCA(obj);
            end
        end

        function obj = LoadData(obj)
            %Loads the data into the app.
            %   When the Excel template is selected, MATLAB reads that file
            %   as a 'table'. Here, the data is separated to be read as
            %   different data types (double, cell, string, etc.). 
            
            tData = readtable(obj.FullPath);
            cVar = tData.Properties.VariableNames;

            %Import again, but make sure all elemental data is numeric
            opts = detectImportOptions(obj.FullPath);
            opts = setvartype(opts, cVar(9:end), 'double');
            tData = readtable(obj.FullPath, opts); 
            %maybe set dates to datetime here?
            tData.SampleDate = datetime(tData.SampleDate);

            % Separate out the data. 
            %Separate out the mixture data and endmember data into MixTable
            %and EMTable
            vDataType = string(tData.DataType);

            %Separate out the mixture data.
            tMixtureData = tData(vDataType == "Mixture", :);

            %Separate out the endmember data.
            tEndmemberData = tData(vDataType == "Endmember", :);

            %Call CleanData from here after raw data is loaded. 
            obj.MixTable = obj.CleanData(tMixtureData);
            obj.EMTable = obj.CleanData(tEndmemberData);
            obj.DataAll = [obj.MixTable; obj.EMTable];
        end

        function obj = FindCat(obj)
            %Separate out each endmember type and average each respective
            %endmember type together. 
            obj.Locations = unique(obj.MixTable.SampleLocation);
            obj.TypeEM = unique(obj.EMTable.SampleLocation);
            obj.NumLoc = length(unique(obj.MixTable.SampleLocation)); 
            obj.NumEM = length(unique(obj.EMTable.SampleLocation));

            %Separate out earliest and latest sample date for all samples. 
            vDatesAll = obj.MixTable.SampleDate;

            obj.StartDate = vDatesAll(1);
            obj.EndDate = vDatesAll(end);

            %Default selected EM. 
            obj.SelectedEM = true(1, length(obj.TypeEM));

            %Default selected location(s).
            obj.SelectedLoc = true(1, length(obj.Locations));

            %Calculate the correlation coefficients for all tracers for
            %just the mixture dataset. 
            obj.CorrMatrix = corrcoef(table2array(obj.MixTable(:, 12:end)));

            %Separate out the names of the array of tracers in dataset as a
            %cell array.
            obj.Tracers = obj.DataAll.Properties.VariableNames(12:end);

            %Default selected tracers. The default is that only the first
            %tracer is selected. This makes generating the correlation
            %matrix in the app much quicker. 
            vSelectedTracers = false(1, length(obj.Tracers));
            vSelectedTracers(1) = 1;
            vSelectedTracers(2) = 1;
            obj.SelectedTracers = vSelectedTracers;
        end

        function obj = SeparateData(obj)
            %Separate out mixture data by location. 
            cMixData = cell(1, obj.NumLoc);

            for i = 1:obj.NumLoc
                cLocation = string(obj.Locations(i));
                mLocData = obj.MixTable(obj.MixTable.SampleLocation == cLocation, :);

                cMixData{i} = mLocData;
            end
            obj.MixtureLoc = cMixData;

            %Separate out endmember data by endmember. 
            cEMData = cell(1, obj.NumEM);

            for i = 1:obj.NumEM
                cEndmember = string(obj.TypeEM(i));
                mEMData = obj.EMTable(obj.EMTable.SampleLocation == cEndmember, :);

                cEMData{i} = mEMData;
            end
            obj.Endmembers = cEMData;  

            %Separate mixture data based on selected location and selected
            %start and end dates. 
            obj.NumSelectLoc = obj.MixtureLoc(obj.SelectedLoc);
            tSelectLoc = obj.NumSelectLoc{1};

            for i = 2:length(obj.NumSelectLoc)
                tSelectLoc = [tSelectLoc; obj.NumSelectLoc{i}];
            end
            
            % % % % % % % Location Selector % % % % % % % %
             
            cLocations = obj.Locations(obj.SelectedLoc);

            for i = 1:length(cLocations)
                 % Create the data inputs for the tables.
                cStartDate = obj.StartDate;
                cEndDate = obj.EndDate;
                
                mDataVal = obj.DataAll;
                calcVal = obj.DataAll.SampleLocation; 
                
                % Sort out the data for the location items 
                val = ismember(calcVal, cLocations{i});
                
                % Sort out the data for the date ranges selected
                dateList = obj.DataAll.SampleDate;
                
                valDate = find(ismember(dateList, cStartDate));
                valDate2 = find(ismember(dateList, cEndDate));
                
                datesYeet = dateList(valDate(1):valDate2(end));
                dateRange = ismember(dateList, datesYeet);
                
                % Compare the location logical and date logical to get another
                % logical that includes both of these variables selected.
                compareData = dateRange & val;
                dataFinal = mDataVal(compareData, :);
                
                cTables{i} = dataFinal;
            end
            
            obj.SelectedLocData = cTables;

            obj.NumSelectLoc = obj.MixtureLoc(obj.SelectedLoc);
            tLocTime = obj.SelectedLocData{1};

            for i = 2:length(obj.NumSelectLoc)
                tLocTime = [tLocTime; obj.SelectedLocData{i}];
            end
            
            obj.SelectedMixDataAll = tLocTime;

            %Get the tracer data.
            mSelectLoc = table2array(tLocTime(:, 12:end));
            mSelectLoc = mSelectLoc(:, obj.SelectedTracers);
            obj.FullMixNumeric = mSelectLoc;

            %Separate out the mixture data by location.
            iRowOn = 0;
            obj.FullNumericLoc = cell(1, length(obj.SelectedLocData));  %% YOU NEED TO CHANGE THIS !!!!!

            for i = 1:length(obj.SelectedLocData)
                iRows = size(obj.SelectedLocData{i}, 1);
                obj.FullNumericLoc{i} = mSelectLoc(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end

            [obj.FullMixSelectedZ, obj.mu, obj.sigma] = zscore(mSelectLoc);

            % Separator for selected endmembers
            cSelectEM = obj.Endmembers(obj.SelectedEM);
            tSelectEM = cSelectEM{1};

            for i = 2:length(cSelectEM)
                tSelectEM = [tSelectEM; cSelectEM{i}];
            end
            
        end

        function obj = CalcPCA(obj)
            %Performs PCA on selected locations and uses generated
            %coefficient matrix to calculate PC scores for selected
            %endmembers. 

            %%%%% Calculate PC scores for the selected locations %%%%%
            % Calculate PCA on full matrix 
            [coeff,PC,latent,tsquared,explained,mu2] = pca(obj.FullMixSelectedZ); 

            %This is important to include so each column is properly
            %annotated. 
            PC = array2table(PC);

            % Build struct of PCA output
            obj.MixPCAData = struct('coeff', coeff, 'score', PC, 'latent', latent, ...
                'tsquared', tsquared, 'explained', explained, 'mu', mu2);

            % Second loop will take the scores calculated from the PCA and
            % store them into individual cells by location 
            iRowOn = 0;
            obj.MixturePCA = cell(1, length(obj.SelectedLocData)); 

            for i = 1:length(obj.SelectedLocData)
                iRows = size(obj.SelectedLocData{i}, 1);
                obj.MixturePCA{i} = PC(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end

            %%%%% Calculate PC scores for the selected endmembers %%%%%
            % Separate out endmember data to include only those endmembers
            % selected.
            cSelectEM = obj.Endmembers(obj.SelectedEM);
            obj.SelectedEMData = cSelectEM;
            tSelectEM = cSelectEM{1};

            for i = 2:length(cSelectEM)
                tSelectEM = [tSelectEM; cSelectEM{i}];
            end

            %Get the tracer data.
            mSelectEM = table2array(tSelectEM(:, 12:end));
            mSelectEM = mSelectEM(:, obj.SelectedTracers);

            %Store numerical endmember data into cell by endmember type
            iRowOn = 0;
            cEMData = cell(1, length(cSelectEM));

            for i = 1:length(cSelectEM)
                iRows = size(cSelectEM{i}, 1);
                cEMData{i} = mSelectEM(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end
            
            obj.SelectedEMNum = cEMData;
            
            %Z-score endmember matrix using same mean and standard deviation
            %as the mixture data.
            iTracers = sum(obj.SelectedTracers, 'all');
            mZEndmember = zeros(size(mSelectEM));
            
            for i = 1:iTracers
                mZEndmember(:, i) = mSelectEM(:, i) - obj.mu(i);
                mZEndmember(:, i) = mZEndmember(:, i)/obj.sigma(i);
            end

            %Perform PCA on z-scored endmember data.
            mEndmemberScore = mZEndmember*coeff;
            PC = mEndmemberScore;
            PC = array2table(PC);
            obj.EMScoreMatrix = PC;

            %Store PC score info into cell array.
            iRowOn = 0;
            obj.EMPCA = cell(1, length(cSelectEM));

            for i = 1:length(cSelectEM)
                iRows = size(cSelectEM{i}, 1);
                obj.EMPCA{i} = PC(iRowOn+1:iRowOn+iRows, :);
                iRowOn = iRowOn + iRows;
            end
            
            % Average the endmembers together to get averaged values. 
            mAvgEM = zeros(length(cEMData), iTracers);
            
            for i = 1:length(cEMData)
                mAvgEM(i, :) = mean(cEMData{i});
            end
            
            %Z-score the averaged endmembers
            obj.AvgEMSelectedZ = zeros(size(mAvgEM));

            for i = 1:iTracers
                obj.AvgEMSelectedZ(:, i) = mAvgEM(:, i) - obj.mu(i);
                obj.AvgEMSelectedZ(:, i) = obj.AvgEMSelectedZ(:, i) / obj.sigma(i);
            end

            %Perform PCA on averaged z-scored EMs. 
            mAvgEMScore = obj.AvgEMSelectedZ*coeff;
            PC = mAvgEMScore;
            PC = array2table(PC);
            obj.AvgEMScore = PC;

            %Store PC score info for averaged endmembers into cell array
            obj.AvgEMPCA = cell(1, length(cSelectEM));

            for i = 1:length(cSelectEM)
                obj.AvgEMPCA{i} = PC(i, :);
            end
        end
    end

    methods (Static)

        function tDataClean = CleanData(tData)
            %Clean the data. 
            %Unique to our dataset, we recieved many concentrations of 0
            %mg/L that still appear as a double but is an unrealistic
            %result, so each of these values are set to NAN to be treated
            %like the other samples that were below the detection limit
            %(i.e. all samples that appeared with a "<" sign in the ICP-MS
            %output data). 

            %Clean mixture data.
            mElData = tData{:, 12:end};
            mElData(mElData == 0) = nan; 

            %Get the minimum value of each column and divide by 2. This
            %will set the values to a reasonably small value without
            %removing the data from the dataset. 
            vMins = min(mElData)/2; 
            vMins = repmat(vMins, size(mElData, 1), 1); 
            mElData(isnan(mElData)) = vMins(isnan(mElData));

            %Put data back into table
            tData{:, 12:end} = mElData;
            tDataClean = tData; 
        end
    end
end