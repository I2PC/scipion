% MIRT_CHECK  check the parameter settings

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
%
% This file is part of the Medical Image Registration Toolbox (MIRT).

function [main,optim]=mirt_check(main,optim,n)

% Default parameters
defmain.similarity='ssd';  % similarity measure (SSD,SAD,CC,RC,MS,CD2,MI)
defmain.okno=16;           % mesh window size
defmain.subdivide = 3;     % number of hierarchical levels
defmain.lambda = 0.01;     % regularization weight, 0 for none
defmain.single=0;          % show the mesh deformation at each iteration
defmain.alpha=0.1;         % similarity measure parameter (e.g. for RC, MS)
defmain.ro=0.9;            % a parameter of MS similarity measure (the assumed correlation)
defmain.MIbins=64;         % Number of bins for the MI similarity measure


defoptim.maxsteps = 300;   % Maximum number of iterations
defoptim.fundif = 1e-6;    % Function tolerance stopping condition
defoptim.gamma = 1;        % Initial step size
defoptim.anneal=0.9;       % Annealing rate


% Check the input options and set the defaults
if n<3,
    main.similarity=defmain.similarity;  % similarity measure
    main.okno=defmain.okno;           % mesh window size
    main.subdivide = defmain.subdivide;     % number of hierarchical levels
    main.lambda = defmain.lambda;     % regularization weight, 0 for none
    main.single=defmain.single;          % show the mesh deformation at each iteration
    main.alpha=defmain.alpha;         % similarity measure parameter
    main.ro=defmain.ro;            % a parameter of MS similarity measure (the assumed correlation)
    main.MIbins=defmain.MIbins;
end
if n<4,
    optim.maxsteps = defoptim.maxsteps;
    optim.fundif = defoptim.fundifdef;
    optim.gamma = defoptim.gamma;
    optim.anneal=defoptim.anneal;
end

if ~isfield(main,'similarity') || isempty(main.similarity), main.similarity = defmain.similarity; end;
if ~isfield(main,'okno') || isempty(main.okno), main.okno = defmain.okno; end;
if ~isfield(main,'subdivide') || isempty(main.subdivide), main.subdivide = defmain.subdivide; end;
if ~isfield(main,'lambda') || isempty(main.lambda), main.lambda = defmain.lambda; end;
if ~isfield(main,'single') || isempty(main.single), main.single = defmain.single; end;
if ~isfield(main,'alpha') || isempty(main.alpha), main.alpha = defmain.alpha; end;
if ~isfield(main,'ro') || isempty(main.ro), main.ro = defmain.ro; end;
if ~isfield(main,'MIbins') || isempty(main.MIbins), main.MIbins = defmain.MIbins; end;


if ~isfield(optim,'maxsteps') || isempty(optim.maxsteps), optim.maxsteps = defoptim.maxsteps; end;
if ~isfield(optim,'fundif') || isempty(optim.fundif), optim.fundif = defoptim.fundifdef; end;
if ~isfield(optim,'gamma') || isempty(optim.gamma), optim.gamma = defoptim.gamma; end;
if ~isfield(optim,'anneal') || isempty(optim.anneal), optim.anneal = defoptim.anneal; end;

% some groupwise options
if ~isfield(main,'group') || isempty(main.group), main.group = 1; end;
if ~isfield(optim,'imfundif') || isempty(optim.imfundif), optim.imfundif = 0.1; end;
if ~isfield(optim,'maxcycle') || isempty(optim.maxcycle), optim.maxcycle = 40; end;
main.cycle=0; main.volume=0; 



