function [uFine, MGParam, errValsFine, kStopFine] = fullMultiGridSolver(MGParam, uExact)
%
% FULLMULTIGRIDSOLVER  FMG driver that starts at coarsest level (0)
% and goes up to the finest (level L).
%
% Outputs:
%   uFine       = final solution on the finest grid (level L)
%   MGParam     = updated struct, including appended levelHistory/colorHistory
%   errValsFine = error/residual from final V-cycle
%   kStopFine   = # of iterations used in that final V-cycle
%

  L       = MGParam.L;       
  nL      = MGParam.nL;      
  xLower  = MGParam.xLower;
  xUpper  = MGParam.xUpper;
  px      = xUpper - xLower;
  hFinest = (xUpper - xLower)/nL;

  %-------------------------------------
  % Construct arrays of n,h per level
  %-------------------------------------
  nLevels = L + 1;  
  nArray  = zeros(nLevels,1);
  hArray  = zeros(nLevels,1);

  nArray(L+1) = nL;          
  hArray(L+1) = hFinest;
  for ell = L-1:-1:0
    nArray(ell+1) = nArray(ell+2)/2;
    hArray(ell+1) = 2*hArray(ell+2);
  end

  fprintf('\n=== Starting Full Multigrid Solver ===\n\n');

  %-------------------------------------
  % Solve on coarsest grid (level 0)
  %-------------------------------------

  fprintf('Solving on Coarsest Level (Level 0)...\n');
  n0 = nArray(1);
  h0 = hArray(1);

  uC = zeros(n0+2,1);
  fC = zeros(n0,1);

  xCentersC = linspace(xLower + h0/2, xUpper - h0/2, n0);
  for i = 1:n0
    xx = xCentersC(i);
    uC(i+1) = sin(2.0*pi*xx/px)^3;
  end
  uC = applyBCs(uC);
  fC = FDOperator(uC,h0);

  MGParam_coarse      = MGParam;
  MGParam_coarse.L    = 0;        
  MGParam_coarse.kMax = 10;       
  [uCsol, MGParam_coarse, ~, ~] = multiGridSolver(uC, fC, h0, MGParam_coarse, uC);

  %-------------------------------------
  % Prolong up to each finer level
  %-------------------------------------
  uPrev    = uCsol;
  MGParamC = MGParam_coarse;  

  for ell = 1:L

    fprintf('\nProlonging from Level %d to Level %d...\n', ell-1, ell);

    nCoarse = nArray(ell);
    nFine   = nArray(ell+1);
    hFine   = hArray(ell+1);

    % Prolong solution from (ell-1)->ell
    uFine = zeros(nFine+2,1);
    uFine = applyBCs(uFine);
    uFine(2:nFine+1) = prolongation(uPrev(2:nCoarse+1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% NEW: Mark this outer step as a "red" transition
    %%% We just jumped from big-grid level (ell-1) up to ell
    MGParamC.levelHistory(end+1) = ell;    % "We've arrived at 'ell' now"
    MGParamC.colorHistory(end+1) = 1;      % 1 => red
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Restrict the exact solution to the current level
    exactFine = uExact;  % Start with the finest-level exact solution
    for restrictLevel = L:-1:ell+1
    % Extract interior points (remove ghost layers)
      exactInterior = exactFine(2:end-1);

    % Apply the restriction operator
      restrictedInterior = restriction(exactInterior);

    % Add ghost layers back to the restricted solution
      exactFine = zeros(length(restrictedInterior) + 2, 1);
      exactFine(2:end-1) = restrictedInterior;
end

% Apply boundary conditions to the restricted solution
exactFine = applyBCs(exactFine);


    % Build the forcing term on the current grid
    fFine = FDOperator(exactFine, hFine);


    % Now do a standard V-cycle on this new level
    MGParam_here   = MGParamC;
    MGParam_here.L = ell;
    [uSol, MGParam_here, errVals, kStop] = multiGridSolver( ...
                                     uFine, fFine, hFine, MGParam_here, exactFine );

    fprintf('Finished V-cycles on Level %d.\n', ell);

    % Save solution & updated MGParam for next step
    uPrev   = uSol;
    MGParamC = MGParam_here;
  end

  %-------------------------------------
  % Return final results at level L
  %-------------------------------------
  uFine       = uPrev;
  MGParam     = MGParamC;
  errValsFine = errVals;
  kStopFine   = kStop;

  fprintf('\n=== Finished Full Multigrid Solver ===\n\n');

end
