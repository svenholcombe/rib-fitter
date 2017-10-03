classdef (Abstract) RibParameterFitterNorm < handle
    
    properties
        targetPointCount = 101; % The number of points to place in the target rib
        ModelXY = nan(1,2);
        TargetXY = nan(1,2);
        PenaltySSEFactors
    end
    properties (Abstract)
        Params;
        Primitives;
    end
    properties (Dependent)
        FreeParams
        BasePrimitives
        NormedModelXY
        NormedTargetXY
    end
    
    methods (Static)
        P = Coefficients_Holcombe2016_JBiomech(this)
    end
    
    
    methods
        function X = get.NormedModelXY(this)
            X = this.ModelXY ./ this.TargetXY(end,1);
        end
        function X = get.NormedTargetXY(this)
            X = this.TargetXY ./ this.TargetXY(end,1);
        end
        function X = get.FreeParams(this)
            X = this.Params(strcmp('free',this.Params.Type),:);
        end
        function set.FreeParams(this,X)
            this.Params(strcmp('free',this.Params.Type),:) = X;
        end
        function X = get.BasePrimitives(this)
            X = this.Primitives(strcmp('base',this.Primitives.Type),:);
        end
        function this = RibParameterFitterNorm
            % When initialising we will set VALUES to NaN
            this.Params.Value(:,1) = nan;
            this.Params.Properties.RowNames = this.Params.Name;
            this.Primitives.XY(:,1) = {zeros(0,2)};
            this.Primitives.Properties.RowNames = this.Primitives.Name;
        end
        function varargout = ParamValues(this,varargin)
            if nargin==2 && iscellstr(varargin{1})
                args = varargin{1};
            else
                args = varargin;
            end
            arrayOut = this.Params.Value(args);
            if nargout>1
                varargout = num2cell(arrayOut);
            else
                varargout = {arrayOut};
            end
        end
    end
    
    methods (Abstract)
        %% IN-PLANE model contraint functions
        [UB,LB] = constraints_ModelParamBounds(this)
        [Aineq, Bineq, Aeq, Beq] = constraints_ModelParamsLinearConstraints(this)
        [cineq, ceq, gineq, geq] = constraints_ModelParamsNonLinearConstraints(ParamsArray,Sx)
        typicalX = initialisation_typicalParameters(Sx)
    end
    methods (Abstract)
        setInherentParameterValues(this);
        setDeterministicParameterValues(this);
        setModelSegmentCoordinates(this,ptSpacing);
        %calculatePenaltySSE(this)
    end
    methods (Static)
        %% IN-PLANE utilities
        function modelXY = resampleAndCleanXYpts(modelXY, ptSpacing)
            % Worker utility to simply take an N-by-2 array of points,
            % remove any NaNs, and resample at ptSpacing equidistant pts.
            modelXY(any(isnan(modelXY),2),:) = [];
            CS = cat(1, 0, cumsum(sqrt(sum(diff(modelXY,[],1).^2,2))));
            repeatMask = cat(1,false,diff(CS)==0);
            modelXY = modelXY(~repeatMask,:);
            CS = CS(~repeatMask);
            modelXY = interp1(CS, modelXY, linspace(0,CS(end), ceil(CS(end)/ptSpacing)));
        end
        function [SSE, minPtDists] = calculate_inPlaneSSE(modelXY, targetXY)
            % What is the in-plane distance between a targetXY and the
            % nearest points on a given modelXY.
            %
            % modelXY and targetXY must be "clean" (no NaNs) and should be
            % close to equi-spaced set of XY points. modelXY should ideally
            % have many more XY points than targetXY.
            %
            % Comparison will occur between each target point and the
            % nearest-neighbour from a window of model points.
            
            modelPtCount = size(modelXY,1);
            targetPtCount = size(targetXY,1);
            modelSearchWindow = ceil(0.5*(modelPtCount/targetPtCount));
            modelSearchCentre = round(linspace(1,modelPtCount,targetPtCount))'; %%%%%%%%%%% FIXME
            searchFrom = min(modelPtCount,max(1,modelSearchCentre-modelSearchWindow));
            searchTo = min(modelPtCount,max(1,modelSearchCentre+modelSearchWindow));

            biggestWindow = max(searchTo - searchFrom);
            compPts = nan(biggestWindow,2,targetPtCount);
            for i = 1:targetPtCount
                compPts(1 + (0:searchTo(i)-searchFrom(i)),:,i) = modelXY(searchFrom(i):searchTo(i),:);
            end
            minSqDists = permute(nanmin(sum(bsxfun(@minus, compPts, permute(targetXY,[3 2 1])).^2,2),[],1),[3 2 1]);
            
            minPtDists = minSqDists; % This is actually sum of SQUARED error
            SSE = sum(minPtDists);
        end
    end
    
    methods
        %% IN-PLANE model fitting
        function inPlaneXY = setAndBuildModelFromParams(this, inParamTable)
            % Set each required param
            paramsToSet = this.Params.Name(ismember(this.Params.Type,{'inherent','free'}));
            for i = 1:length(paramsToSet)
                pName = paramsToSet(i);
                vMatch = strcmpi(pName, inParamTable.Properties.VariableNames);
                vMatchStr = inParamTable.Properties.VariableNames{vMatch};
                this.Params.Value(pName) = inParamTable.(vMatchStr);
            end
            this.setDeterministicParameterValues();
            this.setModelSegmentCoordinates();
            inPlaneXY = this.ModelXY;
        end
        function setAllParameterValues(this,inParams)
            this.setInherentParameterValues;
            this.setOptParameterValues(inParams);
        end
        function setOptParameterValues(this,inParams)
            freeInds = find(strcmp('free',this.Params.Type));
            for i = 1:length(freeInds(:)')
                this.Params.Value(freeInds(i)) = inParams(i);
            end
            this.setDeterministicParameterValues();
        end
        function [thisParams, SSE, flg, History] = getInPlaneOptimisedModelFit(this, init_Params, XYtarget, outputFunctionHandle)
            % Fit [Sh Sz R xR yR] parameters to given N-by-2 XY data and Sx
            % [ShSzRRxRy_init, SSE, flg, stop] = this.getInPlaneOptimisedModelFit(Params0, XY, Sx)
            % ShSzRRxRy_init order = [spSh spSz R xR yR]
            % ShSzRRxRy_init (initial guess) may be empty []
            
            % Register the target
            this.TargetXY = XYtarget;
            this.setInherentParameterValues;
            Sx = this.Params.Value('Sx');
            
            % Apply linear constraints and bounds on parameters
            [Aineq, Bineq, Aeq, Beq] = this.constraints_ModelParamsLinearConstraints();
            [UB,LB] = this.constraints_ModelParamBounds();
            
            % Set up a set of initial guesses based on typical parameters.
            if isempty(init_Params)
                typicalX = this.initialisation_typicalParameters();
            else
                typicalX = init_Params;
            end
            typicalSSE = inf(size(typicalX(:,1)));
            for i = 1:size(typicalX,1)
                thisX = typicalX(i,:);
                constraintsBad = false;
                if any(this.constraints_ModelParamsNonLinearConstraints(thisX)>0) ...
                        || any(thisX>UB | thisX<LB)
                    constraintsBad = true;
                end
                if ~isempty(Aineq) && any(thisX * Aineq' > Bineq)
                    constraintsBad = true;
                end
                % Now we try to get an SSE
                if ~constraintsBad
                    [typicalSSE(i), ~] = this.getInPlaneModelSSE_nonDecimated(thisX, false);
                end
            end
            % Sort the initial guesses by fit
            [~,bestI] = sort(typicalSSE);
            
            % Set up optimisation options
            Opts = optimoptions(@fmincon,'Display','iter', ...
                'Algorithm','interior-point', 'MaxFunEvals', 6000, 'TolX', 1e-8, 'TolFun', 1e-7);
            if nargin >=4 && isa(outputFunctionHandle,'function_handle')
                Opts.OutputFcn = outputFunctionHandle;
            else
                Opts.OutputFcn = @outfun;
            end
            
            % Start the optimisation at the top few initial conditions,
            % breaking when a desirable result is found
            maxOptIterations = 12;
            SSEthresh_exitNow = Sx / 2.5;
            SSEthresh_allEqual = Sx / 2;
            initConCount = min(maxOptIterations,size(typicalX,1));
            outParams = nan(initConCount, size(typicalX,2));
            SSEs = nan(initConCount,1); Histories = cell(initConCount,1);
            
            for i = 1:initConCount
                History = struct; % Initialise empty history output
                init_Params = typicalX(bestI(i),:);
                [thisParams, SSE, flg] = fmincon(@(currParams)this.getInPlaneModelSSE_nonDecimated(currParams, true),...
                    init_Params, Aineq, Bineq(:), Aeq, Beq, LB, UB, @this.constraints_ModelParamsNonLinearConstraints, Opts);
                outParams(i,:) = thisParams; SSEs(i) = SSE; Histories{i} = History;
                % Break under 3 conditions: firstly if a really low SSE was
                % found. Secondly if multiple OK optimisations all reached
                % the same result. Thirdly if we've tried 5 times.
                if SSE < SSEthresh_exitNow
                    break;
                elseif i>2 && SSE<SSEthresh_allEqual && ~any(abs(diff(SSEs(1:i)))>1.0)
                    break;
                end
            end
            [~,minIdx] = min(SSEs);
            History = Histories{minIdx};
            thisParams = outParams(minIdx,:);
            % Run one last time using the best params to make sure that the
            % fitter object's this.Params match the best result;
            SSE = this.getInPlaneModelSSE_nonDecimated(thisParams, true);
            this.setDeterministicParameterValues;
            this.setOutcomeParameters;
            this.setModelSegmentCoordinates(); % Unneeded unless wanting to use this.plot()
            
            function stop = outfun(x,optimValues,state)
                stop = false;
                switch state
                    case 'init'
                        History = struct('x',x,'fval',optimValues.fval);
                    case 'iter'
                        % Concatenate current point and objective function
                        % value with history. x must be a row vector.
                        History.fval = [History.fval; optimValues.fval];
                        History.x = [History.x; x];
                    otherwise
                end
            end
        end
    end
end