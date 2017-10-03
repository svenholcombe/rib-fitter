classdef RibModeler < RibParameterFitterNorm
    properties
        Params = cell2table({
            'Sx'            'S_x'           'inherent'    'mm'      'Rib Length'
            'PkX'           'X_{Pk}'        'free'        'mm/S_x'  'Rib peak location (X)'
            'PkY'           'Y_{Pk}'        'free'        'mm/S_x'  'Rib peak location (Y)'
            'proxInnAngle'  '\phi_{pia}'    'free'        'degs'    'Proximal Inner Angle (degs)'
            'Bp'            'B_{p}'         'free'        '.'       'Proximal log-spiral shape'
            'Bd'            'B_{d}'         'free'        '.'       'Distal log-spiral shape'
            'proxThPk'      '\theta{}_p'    'determined'  'rads'    'Theta (rads) at which proximal spiral reaches its end at the rib peak'
            'proxThO'       '\theta{}_o'    'determined'  'rads'    'Theta (rads) at which proximal spiral reaches its end at the origin'
            'proxThCh'      '\theta{}_{ch}' 'determined'  'rads'    'Theta (rads) at which proximal spiral reaches its convex hull to the origin'
            'distThEnd'     '\theta{}_d'    'determined'  'rads'    'Theta (rads) for distal spiral at rib distal end'
            'distThPk'      '\theta{}_x'    'determined'  'rads'    'Theta (rads) for distal spiral at its peak Y coordinate'
            'distInnAngle'  '\phi_{dia}'    'outcome'     'degs'    'Inner Angle (degs) between the distal end of the rib and the x-axis'
            'lenCentroid2d' 'L_{2d}'        'outcome'     'mm'      'Centroidal model length of the rib'
            'distKEnd'      '\kappa{}_{dist}'   'outcome' 'mm^{-1}' 'Curvature (inverse of radius of curvature) at distal end'
            'proxKCh'       'K_{proxCH}'    'outcome'     'mm^{-1}' 'Curvature (inverse of radius of curvature) at proximal convex hull tangent point'
            'proxKPost'     '\kappa{}_{post}'  'outcome'  'mm^{-1}' 'Curvature (inverse of radius of curvature) at proximal end at location of maximum posterior extrusion'
            'fit_spsse'     'SSE'           'fitresult'   'mm^2'    'In-plane error (sum of squared distances to target)'
            'fit_meanerr'   '\mu Err'       'fitresult'   'mm'      'In-plane error (mean of abs distances to target)'
            'ang_ph'        '\alpha_{PH}'   'mdr'         'degs'    'Pump-handle angle'
            'ang_ls'        '\alpha_{LS}'   'mdr'         'degs'    'Lateral swing angle'
            'ang_bh'        '\alpha_{BH}'   'mdr'         'degs'    'Bucket-handle angle'
            },'Var',{'Name','LabelLatex','Type','Unit','Definition'});
        
        Primitives = cell2table({
            'Origin2Prox'   'join'      {}                  [    0     0     0]
            'DistSpiral'    'base'      {'Bd','PkX','PkY'}  [0.494 0.184 0.556]
            'ProxSpiral'    'base'      {'Bp','PkX','PkY'}  [0.466 0.674 0.188]
            },'Var',{'Name','Type','Params','RGB'});
    end
    
    methods
        %% IN-PLANE model contraint functions
        function [UB,LB] = constraints_ModelParamBounds(~)
            % Apply bounds on parameters
            %    [PkX,  PkY, PIA,    Bp,   Bd]
            UB = [0.8   0.9  160,   3.0,  3.5];  %upper bounds
            LB = [0.05 0.05    1,  -3.0, -3.5];  %lower bounds
        end
        function [Aineq, Bineq, Aeq, Beq] = constraints_ModelParamsLinearConstraints(~)
            % Linear constraints to parameters [Bp Bd]
            [Aineq,Bineq] = deal(zeros(0,2));
            Aeq=[]; Beq=[];             %no equality constraints
        end
        function [cineq, ceq] = constraints_ModelParamsNonLinearConstraints(this,inPset)
            % Non-linear constraints to parameters [PkX PkY PIA Bp Bd]
            ceq = [];
            
            PkX = inPset(1); PkY = inPset(2); PIA = inPset(3); Bp = inPset(4);
            %% We will need functions for x, y, and 1st derivatives
            xfcn = @(t) exp(Bp.*t) .* cos(t);
            yfcn = @(t) exp(Bp.*t) .* sin(t);
            dxdtfcn = @(t)Bp*exp(Bp.*t).*cos(t) - exp(Bp.*t).*sin(t);
            dydtfcn = @(t)Bp*exp(Bp.*t).*sin(t) + exp(Bp.*t).*cos(t);
            
            %% The PeakXY location defines an angle
            
            % There will be a theta at which the spiral's gradient matches
            % the required gradient to meet [PkX,PkY] with zero slope
            dydxAtan2 = atan2d(PkY,PkX);
            % Solve for that first cut location (first theta, then get X,Y)
            cutTh = fzero(@(t)atan2d(-dydtfcn(t),-dxdtfcn(t)) + dydxAtan2,pi/2);
            cutx = xfcn(cutTh); %#ok<NASGU>
            cuty = yfcn(cutTh);
            
            % We are looking for the second cut. It will occur on the
            % opposite side to cut1, between the theta at spiral max and
            % max+pi (which is the min of the spiral)
            thetaAtMax = 2 * atan(sqrt(Bp^2+1)+Bp);
            thetaAtMin = thetaAtMax + pi;
            % It will be the location matching the original cut's Y
            % coordinate. However, if that can't be matched because the
            % spiral doesn't go down that far on the second cut side, we'll
            % just use the closest cut possible (BUT WE SHOULD PENALIZE
            % THIS EVENT SO THAT A POOR CUT LOCATION ISN'T CHOSEN!)
            if cuty > yfcn(thetaAtMin)
                oppCutTh = fzero(@(t)yfcn(t)-cuty,[thetaAtMax thetaAtMin]);
                this.PenaltySSEFactors = [];
            else
                oppCutTh = thetaAtMin;
                this.PenaltySSEFactors = 2.0;
            end
            
            angDtangent = PIA - dydxAtan2;
            % The tangential angle at the opposite cut location is:
            angDoppCutTh = atan2d(-dydtfcn(oppCutTh),-dxdtfcn(oppCutTh));
            cineq = [
                yfcn(thetaAtMin) - cuty
                angDtangent - angDoppCutTh
                dydxAtan2 - PIA + 1];
        end
        function typicalX = initialisation_typicalParameters(~)
            %      PkX       PkY       PIA        Bp        Bd
            typicalX = [
                0.3791    0.5241   80.9864   -0.5694    1.6813
                0.3212    0.6199  111.0250   -0.3665    0.7061
                0.3432    0.5885  116.1477   -0.5106    0.3805
                0.3247    0.5295  114.0817   -0.4978   -0.0010
                0.3179    0.4853  112.1643   -0.5142   -0.3090
                0.3074    0.4455  110.1416   -0.4634   -0.4294
                0.2980    0.3993  107.5676   -0.3880   -0.2146
                0.3163    0.3716  103.9950   -0.4000   -0.1657
                0.3599    0.3618   96.0719   -0.5321   -0.0097
                0.3788    0.3347   89.8284   -0.5940   -0.0745
                0.3941    0.2801   84.7475   -0.6546    0.3707
                0.4076    0.1772   55.5660   -1.7505    1.5303
                0.3206    0.7824  103.4729   -0.0656    0.3695
                0.2445    0.3253  135.9036   -0.5177   -0.4853
                0.2536    0.3323  148.0421   -0.6316   -1.6761
                0.2272    0.2000   75.3246   -0.8136    0.2690
                0.2000    0.7000  116.8506         0    0.2690
                0.2000    0.6000  112.6199         0   -0.0400
                0.2200    0.3300  116.1250   -0.6650   -1.3380
                0.2970    0.5386  135.8224   -0.5801   -0.3073
                0.1329    0.5477  152.7280   -0.2427   -0.3179
                ];
        end
    end
    
    properties (Hidden)
        BdInvFcnX = @(X)X;
        BdInvFcnY = @(Y)Y;
        ProxInvFcnXYtoTheta = @(XY)XY;
        DistInvFcnXYtoTheta = @(XY)XY;
        DistFcnTangAngByTheta = @(t)t;
        DistMagRatio = nan;
        ProxFcnTangAngByTheta = @(t)t;
        ProxMagRatio = nan;
        
        
        BdFcnXbyTheta = @(t)t;
        BdFcnYbyTheta = @(t)t;
        BpFcnXYbyTheta = @(t)t;
        BdSpiralPtsXY = zeros(0,2);
        BpSpiralPtsXY = zeros(0,2);
    end
    
    methods
        function setInherentParameterValues(this)
            % Based on this.TargetXY, set any "inherent" params
            Sx = this.TargetXY(end,1);
            this.Params.Value('Sx') = Sx;
        end
        function setDistSpiralDeterminedParameters(this)
            % Set the determined values of the spiral, and set function
            % handles to generate X and Y spiral coordinates.
            %%
            % The following parameters are GIVEN
            [Bd,PkX,PkY] = this.ParamValues({'Bd','PkX','PkY'});
            
            % From these, we can derive a logarithmic spiral that meets the
            % above conditions. The first step is to find the theta
            % parameter at which a log-spiral with this shape reaches its
            % "peak" Y-value. This will allow us to relate this shape to
            % the pre-provided PeakX and PeakY coordinates.
            
            % Derivative:
            % (d)/(dt)(e^(Bd t) sin(t)) = e^(Bd t) (Bd sin(t)+cos(t))
            % Maximum occurs at d/dt = 0:
            thetaMax = 2 * atan(sqrt(Bd^2+1)+Bd);
            %
            % Gather the peak in regular spiral XY space
            xFcn = @(t)-exp(Bd*t).*cos(t);
            yFcn = @(t)exp(Bd*t).*sin(t);
            pX = xFcn(thetaMax);
            pY = yFcn(thetaMax);
            % Calculate the required slope from the peak to the rib end
            pSlope = -PkY / (1-PkX);
            % Calculate the theta corresponding to a point on the spiral
            % with that slope from the spiral maximum
            thetaEnd = fzero(@(t)(exp(Bd*t).*sin(t) - pY) ./ (-exp(Bd*t).*cos(t) - pX) - pSlope, pi);
            % Get the X, Y coordinates of that end (in spiral XY space)
            eX = -exp(Bd*thetaEnd).*cos(thetaEnd);
            eY = exp(Bd*thetaEnd).*sin(thetaEnd);
            % Calculate the rescaling ratio to size the spiral to the
            % normalized space
            pMagnitude = sqrt((pX-eX)^2 + (pY-eY)^2);
            rMagnitude = sqrt(PkY^2 + (1-PkX)^2);
            pMagRatio = rMagnitude / pMagnitude;
            
            % Display
            %eeX = @(t)(-exp(Bd*t).*cos(t)-pX) * pMagRatio + PkX;
            %eeY = @(t)( exp(Bd*t).*sin(t)-pY) * pMagRatio + PkY;
            %thetas = linspace(thetaMax, endThP, 50);
            %figure, plot(eeX(thetas),eeY(thetas),'-',[PkX 1],[PkY 0],'-'), axis image
            
            % Set parameters for storage
            this.Params.Value('distThEnd') = thetaEnd;
            this.Params.Value('distThPk') = thetaMax;
            
            % Set functions for normalized XY coordinate placement
            this.BdFcnXbyTheta = @(t)(-exp(Bd*t).*cos(t)-pX) * pMagRatio + PkX;
            this.BdFcnYbyTheta = @(t)( exp(Bd*t).*sin(t)-pY) * pMagRatio + PkY;
            % Set inverse functions
            this.BdInvFcnX = @(X)(X-PkX)/pMagRatio + pX;
            this.BdInvFcnY = @(Y)(Y-PkY)/pMagRatio + pY;
            % Set inverse functions
            this.DistInvFcnXYtoTheta = @(XY)pi - cart2pol((XY(:,1)-PkX)/pMagRatio + pX,(XY(:,2)-PkY)/pMagRatio + pY);
            this.DistFcnTangAngByTheta = @(t)atan2(...
                Bd*exp( Bd*t).*sin(t) + exp(Bd*t).*cos(t),...
                -Bd*exp(Bd*t).*cos(t) + exp(Bd*t).*sin(t));
            this.DistMagRatio = pMagRatio;
            
            %% Some tests at making sure the inverse fcns work
%             t = linspace(thetaMax, thetaEnd, 20)';
%             x = xFcn(t); y = yFcn(t);
%             X = this.BdFcnXbyTheta(t);
%             Y = this.BdFcnYbyTheta(t);
%             XY = [X Y];
%             invx = (XY(:,1)-PkX)/pMagRatio + pX;
%             invy = (XY(:,2)-PkY)/pMagRatio + pY;
%             invC = cart2pol(invx,invy);
%             invT = this.DistInvFcnXYtoTheta([X Y]);
%             if any(abs([invx invy] - [x y])>0.000001)
%                 
%                 XYt = this.DistInvFcnXYtoTheta([X,Y]);
%                 figure, subplot(1,2,1)
%                 plot(x,y , '.'), title('x,y'), axis image, hold on, plot(invx,invy,'o')
%                 subplot(1,2,2)
%                 plot(X,Y , '.'), title('X,Y'), axis image
%                 uiwait
%             end
        end
        function setProxSpiralDeterminedParameters(this)
            % Set the determined values of the proximal spiral, and set
            % function handles to generate X and Y spiral coordinates.
            
            % The following parameters are GIVEN
            [PIA,Bp,PkX,PkY] = this.ParamValues({'proxInnAngle','Bp','PkX','PkY'});
            %% We will need functions for x, y, and 1st derivatives
            xfcn = @(t) exp(Bp.*t) .* cos(t);
            yfcn = @(t) exp(Bp.*t) .* sin(t);
            dxdtfcn = @(t)Bp*exp(Bp.*t).*cos(t) - exp(Bp.*t).*sin(t);
            dydtfcn = @(t)Bp*exp(Bp.*t).*sin(t) + exp(Bp.*t).*cos(t);
            
            %% The PeakXY location defines an angle
            
            % There will be a theta at which the spiral's gradient matches
            % the required gradient to meet [PkX,PkY] with zero slope
            dydxAtan2 = atan2d(PkY,PkX);
            % Solve for that first cut location (first theta, then get X,Y)
            cutTh = fzero(@(t)atan2d(-dydtfcn(t),-dxdtfcn(t)) + dydxAtan2,pi/2);
            cutx = xfcn(cutTh);
            cuty = yfcn(cutTh);
            
            % We are looking for the second cut. It will occur on the
            % opposite side to cut1, between the theta at spiral max and
            % max+pi (which is the min of the spiral)
            thetaAtMax = 2 * atan(sqrt(Bp^2+1)+Bp);
            thetaAtMin = thetaAtMax + pi;
            % It will be the location matching the original cut's Y
            % coordinate. However, if that can't be matched because the
            % spiral doesn't go down that far on the second cut side, we'll
            % just use the closest cut possible (BUT WE SHOULD PENALIZE
            % THIS EVENT SO THAT A POOR CUT LOCATION ISN'T CHOSEN!)
            if cuty > yfcn(thetaAtMin)
                oppCutTh = fzero(@(t)yfcn(t)-cuty,[thetaAtMax thetaAtMin]);
                this.PenaltySSEFactors = [];
            else
                oppCutTh = thetaAtMin;
                this.PenaltySSEFactors = 2.0;
            end
            
            % Based on PkX, PkY and PIA, we can know the tangent to the
            % spiral at the point where we transition to a straight line.
            % Get the location on the spiral with this tangent.
            % The tangential angle required for this PkX, PkY, PIA is:
            angDtangent = PIA - dydxAtan2;
            % The tangential angle at the opposite cut location is:
            angDoppCutTh = atan2d(-dydtfcn(oppCutTh),-dxdtfcn(oppCutTh));
            % We only have real tangent if angDtang occurs above the cut:
            if angDoppCutTh >= angDtangent && angDtangent>0
                tangTh = fzero(@(t)atan2d(-dydtfcn(t),-dxdtfcn(t)) - angDtangent,[thetaAtMax, oppCutTh]);
            else % The tangent location is past the cut. Infeasible
                tangTh = oppCutTh;
                this.PenaltySSEFactors = 1.0 + (angDtangent-angDoppCutTh);
            end
            
            % Get the position of transition from spiral to tangent
            tangx = xfcn(tangTh); tangy = yfcn(tangTh);
            % We know where this tangent crosses the cut line (y==cuty)
            A = tangy - cuty;
            B = A / tand(angDtangent);
            tangEndx = tangx - B;
            
            % A prox scale factor sizes from dist between ends of spiral +
            % tangent in original "x" space and normalized rib "X" space
            lenInx = cutx - tangEndx;
            lenInX = sqrt(PkX^2 + PkY^2);
            proxSF = lenInX / lenInx;
            
            % Build the transformation to shift/scale/rotate [x y]
            rotAng = atan(PkY/PkX);
            rotXform = [cos(rotAng) sin(rotAng) 0; -sin(rotAng) cos(rotAng) 0; 0 0 1];
            scaXform = [proxSF 0 0; 0 proxSF 0; 0 0 1];
            offXform = [1 0 0; 0 1 0; -tangEndx -cuty 1];
            TX = affine2d(offXform * scaXform * rotXform);
            XYfcn = @(t)TX.transformPointsForward([xfcn(t(:)) yfcn(t(:))]);
    
            % Set the inverse functions to transform BACK from XY to xy
            XYtoxyInvFcn = @(XY)TX.transformPointsInverse(XY);
            cartToPolXYarr = @(xy)cart2pol(xy(:,1),xy(:,2)); % Extra inline needed to send X,Y to cart2pol
            this.ProxInvFcnXYtoTheta = @(XY)cartToPolXYarr(XYtoxyInvFcn(XY));
            
            %thetas = linspace(cutTh, oppCutTh, 100)';
            %XY = XYfcn(thetas);  XYch = XYfcn(chTh);
            %figure, plot(XY(:,1),XY(:,2),'-',[0 XYch(1)],[0 XYch(2)],'.-'), axis image, grid on
            %%
            this.BpFcnXYbyTheta = XYfcn;
            this.Params.Value('proxThPk') = cutTh;
            this.Params.Value('proxThO') = oppCutTh;
            this.Params.Value('proxThCh') = tangTh;
            
            this.ProxFcnTangAngByTheta = @(t)atan2(dydtfcn(t),dxdtfcn(t)) + atan2(PkY,PkX);
            this.ProxMagRatio = proxSF;
        end
        
        %% IN-PLANE model fitting
        function [SSE, minPtDistsSq] = getInPlaneModelSSE_nonDecimated(this, inParams, applyPenalties)
            % First set parameters
            freeInds = find(strcmp('free',this.Params.Type));
            for i = 1:length(freeInds(:)')
                this.Params.Value(freeInds(i)) = inParams(i);
            end
            if any(isnan(inParams))
                SSE = 1e6;
                return;
            end
            % Make sure spiral deterministic params are fully calculated
            this.setProxSpiralDeterminedParameters();
            this.setDistSpiralDeterminedParameters();
            
            % Then measure distances from target to model
            [SSE, minPtDistsSq] = calculate_inPlaneSSEmodelToTarget(this);
            % Apply penalties if requested
            if nargin>=3 && applyPenalties
                for i = 1:length(this.PenaltySSEFactors)
                    SSE = SSE + 1000; %this.PenaltySSEFactors(i) * SSE;%this.calculatePenaltySSE();
                end
            end
        end
        function [SSE, minPtDistsSq] = calculate_inPlaneSSEmodelToTarget(this)
            %%
            [Sx,PkX,proxThCh] = this.ParamValues('Sx','PkX','proxThCh');
            
            ribXYnormed = this.TargetXY ./ Sx;
            distMask = ribXYnormed(:,1) >= PkX;
            
            % Distal side is straight forward
            distXYnormed = ribXYnormed(distMask,:);
            distXYthetas = this.DistInvFcnXYtoTheta(distXYnormed);
            remappedDistXY = [this.BdFcnXbyTheta(distXYthetas) this.BdFcnYbyTheta(distXYthetas)];
            % Proximal side has two portions. Get thetas from all
            proxXYnormed = ribXYnormed(~distMask,:);
            proxXYthetas = this.ProxInvFcnXYtoTheta(proxXYnormed);
            % Now get prox counterparts belonging to spiral
            proxSpiMask = proxXYthetas <= proxThCh;
            remappedProxXY = this.BpFcnXYbyTheta(proxXYthetas(proxSpiMask));
            % And prox counterparts belonging to line
            proxChXY = this.BpFcnXYbyTheta(proxThCh);
            
            % Any points belonging to the line can be distanced directly
            % from a line running from origin to CH point
            proxLineDists = distancePointEdge(proxXYnormed(~proxSpiMask,:),[0 0 proxChXY]);
            
            %distSpiDists = sqrt(sum((distXYnormed - remappedDistXY).^2,2));
            %proxSpiDists = sqrt(sum((proxXYnormed(proxSpiMask,:) - remappedProxXY).^2,2));
            % The pt2pt dists above are better calculated by taking each
            % point on the spiral as a line (at that point with direction
            % equal to the spiral tangent at that point) and calculate the
            % projection distance from point to line. This avoids the minor
            % additional error you get if the projected/remapped points are
            % not basically remapped perpendicular to the spiral.
            distAngs = this.DistFcnTangAngByTheta(distXYthetas);
            distPtLines = [remappedDistXY cos(distAngs) sin(distAngs)];
            distSpiDists = distancePointLine(distXYnormed, distPtLines);
            
            proxAngs = this.ProxFcnTangAngByTheta(proxXYthetas(proxSpiMask));
            proxPtLines = [remappedProxXY cos(proxAngs) sin(proxAngs)];
            proxSpiDists = distancePointLine(proxXYnormed(proxSpiMask,:), proxPtLines);
            
            minPtDists = [proxLineDists; proxSpiDists; distSpiDists]*Sx;
            minPtDistsSq = minPtDists.^2;
            SSE = sum(minPtDistsSq);
            this.Params.Value('fit_spsse') = SSE;
            this.Params.Value('fit_meanerr') = mean(minPtDists);
        end
        function setDeterministicParameterValues(this)
            % Make sure spiral deterministic params are fully calculated
            this.setProxSpiralDeterminedParameters();
            this.setDistSpiralDeterminedParameters();
            % Store a copy of a 2000 pt spiral to save computation
            [Sx,distThEnd,distThPk] = this.ParamValues({'Sx','distThEnd','distThPk'});
            t = linspace(distThPk,distThEnd,200)';
            this.BdSpiralPtsXY = [this.BdFcnXbyTheta(t) this.BdFcnYbyTheta(t)]*Sx;
            %
            [proxThCh,proxThPk] = this.ParamValues({'proxThCh','proxThPk'});
            t = linspace(proxThCh,proxThPk,100)';
            this.BpSpiralPtsXY = this.BpFcnXYbyTheta(t)*Sx;
        end
        function setModelSegmentCoordinates(this)
            %%
            [Sx,proxThCh] = this.ParamValues('Sx','proxThCh');
            minPtSpacing = 0.1 / Sx;
            
            % Domain from the origin to the proximal spiral convex hull pt
            this.Primitives.XY{'Origin2Prox'} = [0 0; this.BpFcnXYbyTheta(proxThCh)]*Sx;
            % Domain from the prox spiral convex hull to the peak
            this.Primitives.XY{'ProxSpiral'} = this.BpSpiralPtsXY;
            % Domain from the peak to the end of the distal spiral
            this.Primitives.XY{'DistSpiral'} = this.BdSpiralPtsXY;
            % Connect all domains together
            modelXY = [
                this.Primitives.XY{'Origin2Prox'}
                this.Primitives.XY{'ProxSpiral'}
                this.Primitives.XY{'DistSpiral'}];
            
            normedModelXYpts = this.resampleAndCleanXYpts(modelXY, minPtSpacing) / Sx;
            CS = cat(1,0,cumsum(sqrt(sum(diff(normedModelXYpts,[],1).^2,2))));
            this.ModelXY = interp1(CS, normedModelXYpts, linspace(0,CS(end), this.targetPointCount)) * Sx;
        end
        
        function setOutcomeParameters(this)
            
            % Fetch some relevant parameters
            [Bd,Sx,distThEnd,Bp,proxThCh] = this.ParamValues('Bd','Sx','distThEnd','Bp','proxThCh');
            
            % The inner angle is the negative of the slope at distal end
            this.Params.Value('distInnAngle') = rad2deg(this.DistFcnTangAngByTheta(distThEnd) * -1);
            % And the inverse of the slope at the proximal end (conv hull)
            this.Params.Value('proxInnAngle') = 180 + rad2deg(this.ProxFcnTangAngByTheta(proxThCh));
            
            % Curvature of a spiral is based on the simple equation:
            kFcnBase = @(a,B,t)exp(-B*t) ./ (a*sqrt(1+B^2));
            % The curvature can be interrogated at the distal end
            this.Params.Value('distKEnd') = kFcnBase(this.DistMagRatio,Bd,distThEnd) / Sx;
            % And the curvature can be interrogated at the proximal end
            this.Params.Value('proxKCh') = kFcnBase(this.ProxMagRatio,Bp,proxThCh) / Sx;
            % The posterior-most location will be the location of local
            % maximum given by the equation:
            thetaAtMax = 2 * atan(sqrt(Bp^2+1)+Bp);
            this.Params.Value('proxKPost') = kFcnBase(this.ProxMagRatio,Bp,thetaAtMax) / Sx;
            
            % We'll just calculate the length via the sum of point dists
            this.Params.Value('lenCentroid2d') = sum(sqrt(sum(diff(this.ModelXY(:,1:2)).^2,2)));
        end
        
        function varargout = getInPlaneOptimisedModelFit(this,init_Params,varargin)
            [thisParams, SSE, flg, History] = getInPlaneOptimisedModelFit@RibParameterFitterNorm(this,init_Params,varargin{:});
            varargout = {thisParams, SSE, flg, History};
            if isempty(init_Params)
                % Make a new set of init params with added 15 deg PIA and
                % drop the Bp a little
                new_init = thisParams;
                new_init(3) = new_init(3) + 15;
                new_init(4) = new_init(4) * 0.5;
                [tp2, SSE2, flg2, History2] = getInPlaneOptimisedModelFit@RibParameterFitterNorm(this,new_init,varargin{:});
                % If this new init conditions are better, use them
                if SSE2 < (SSE - 1e-5)
                    fprintf('***********\nBetter fit found from -15 deg refit!\n')
                    fprintf('Old: [%s] (SSE=%0.1f)\nNew: [%s] (SSE=%0.1f)\n***********\n',...
                        num2str(thisParams), SSE, num2str(tp2), SSE2)
                    varargout = {tp2, SSE2, flg2, History2};
                end
            end
            % Run one last time using the best params to make sure that the
            % fitter object's this.Params match the best result;
            varargout{2} = this.getInPlaneModelSSE_nonDecimated(varargout{1}, true);
        end
        
        function plot(this, varargin)

            IP = inputParser;
            IP.addParameter('Parent',[]);
            IP.addParameter('Normed',false);
            IP.parse(varargin{:})
            axH = IP.Results.Parent;
            if isempty(axH) || ~ishandle(axH)
                axH = gca;
            end
            plot2 = @(x,varargin)plot(x(:,1),x(:,2),varargin{:});
            [Bp,PIA,Bd,PkX,PkY,Sx] = this.ParamValues({'Bp','proxInnAngle','Bd','PkX','PkY','Sx'});
            
            % Make a scale factor for unnormed coords (SFU) and normed (SF)
            SF = Sx;
            if IP.Results.Normed
                SF = 1;
            end
            SFU = SF / Sx;
            mXY = this.ModelXY * SFU;
            mSegsXY = cellfun(@(x)x*SFU, this.Primitives.XY,'Un',0);
            
            pXY = this.BpSpiralPtsXY * SFU;
            dXY = this.BdSpiralPtsXY * SFU;
            hold(axH,'on')
            plot2(mXY,'LineWidth',2)
            pSpiH = plot2(pXY,'--','Color',this.Primitives.RGB('ProxSpiral',:));
            plot2(dXY,'--','Color',this.Primitives.RGB('DistSpiral',:));
            plot2(this.Primitives.XY{'Origin2Prox'}*SFU,'.k-')
            
            plot(PkX*SF,PkY*SF,'k.')
            cellfun(@(x)plot2(x([1 end],:),'k.'), mSegsXY(~cellfun(@isempty,mSegsXY)))
            axis equal
            text(axH.XLim(2),axH.YLim(1),sprintf('(Bp,PIA) = (%0.3f, %0.1f), Bd = %0.3f, PkXY = [%0.2f, %0.2f]', ...
                [Bp PIA, Bd PkX PkY]),'Hor','Right','Vert','bottom',...
                'Color',pSpiH.Color)
        end
    end
end