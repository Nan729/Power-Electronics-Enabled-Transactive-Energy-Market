function mpc = case33open_new_DG

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
    1   3   0       0       0   0   1   1.05   0   10   1   1.1   0.9;
	2   1   0.100   0.060   0   0   1   1   0   10   1   1.1   0.9;
    3   1   0.090   0.040   0   0   1   1   0   10   1   1.1   0.9;
    4   1   0.120   0.080   0   0   1   1   0   10   1   1.1   0.9;
    5   1   0.060   0.030   0   0   1   1   0   10   1   1.1   0.9;
    6   1   0.060   0.020   0   0   1   1   0   10   1   1.1   0.9;
    7   1   0.200   0.100   0   0   1   1   0   10   1   1.1   0.9;
    8   1   0.200   0.100   0   0   1   1   0   10   1   1.1   0.9;
    9   1   0.060   0.020   0   0   1   1   0   10   1   1.1   0.9;
   10   1   0.060   0.020   0   0   1   1   0   10   1   1.1   0.9;
   11   1   0.045   0.030   0   0   1   1   0   10   1   1.1   0.9;
   12   1   0.060   0.035   0   0   1   1   0   10   1   1.1   0.9;
   13   1   0.060   0.035   0   0   1   1   0   10   1   1.1   0.9;
   14   1   0.120   0.080   0   0   1   1   0   10   1   1.1   0.9;
   15   1   0.060   0.010   0   0   1   1   0   10   1   1.1   0.9;
   16   1   0.060   0.020   0   0   1   1   0   10   1   1.1   0.9;
   17   1   0.060   0.020   0   0   1   1   0   10   1   1.1   0.9;
   18   1   0.090   0.040   0   0   1   1   0   10   1   1.1   0.9;
   19   1   0.090   0.040   0   0   1   1   0   10   1   1.1   0.9;
   20   1   0.090   0.040   0   0   1   1   0   10   1   1.1   0.9;
   21   1   0.090   0.040   0   0   1   1   0   10   1   1.1   0.9;
   22   1   0.090   0.040   0   0   1   1   0   10   1   1.1   0.9;
   23   1   0.090   0.050   0   0   1   1   0   10   1   1.1   0.9;
   24   1   0.420   0.200   0   0   1   1   0   10   1   1.1   0.9;
   25   1   0.420   0.200   0   0   1   1   0   10   1   1.1   0.9;
   26   1   0.060   0.025   0   0   1   1   0   10   1   1.1   0.9;
   27   1   0.060   0.025   0   0   1   1   0   10   1   1.1   0.9;
   28   1   0.060   0.020   0   0   1   1   0   10   1   1.1   0.9;
   29   1   0.120   0.070   0   0   1   1   0   10   1   1.1   0.9;
   30   1   0.200   0.600   0   0   1   1   0   10   1   1.1   0.9;
   31   1   0.150   0.070   0   0   1   1   0   10   1   1.1   0.9;
   32   1   0.210   0.100   0   0   1   1   0   10   1   1.1   0.9;
   33   1   0.060   0.040   0   0   1   1   0   10   1   1.1   0.9;
];

%% branch data
%	fbus	tbus	r(Ohm)	x(Ohm)	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    1   2   0.0922   0.0470   0   250   250   250   0   0   1;
    2   3   0.4930   0.2511   0   250   250   250   0   0   1;
    3   4   0.3660   0.1864   0   250   250   250   0   0   1;
    4   5   0.3811   0.1941   0   250   250   250   0   0   1;
    5   6   0.8190   0.7070   0   250   250   250   0   0   1;
    6   7   0.1872   0.6188   0   250   250   250   0   0   1;
    7   8   0.7114   0.2351   0   250   250   250   0   0   1;
    8   9   0.3300   0.7400   0   250   250   250   0   0   1;
    9   10   0.2440   0.7400   0   250   250   250   0   0   1;
   10   11   0.1966   0.0650   0   250   250   250   0   0   1;
   11   12   0.3744   0.1238   0   250   250   250   0   0   1;
   12   13   0.4680   1.1550   0   250   250   250   0   0   1;
   13   14   0.5416   0.7129   0   250   250   250   0   0   1;
   14   15   0.5910   0.5260   0   250   250   250   0   0   1;
   15   16   0.6463   0.5450   0   250   250   250   0   0   1;
   16   17   0.5890   1.7210   0   250   250   250   0   0   1;
   17   18   0.6320   0.5740   0   250   250   250   0   0   1;
    2   19   0.1640   0.1565   0   250   250   250   0   0   1;
   19   20   0.6042   1.3554   0   250   250   250   0   0   1;
   20   21   0.4095   0.4784   0   250   250   250   0   0   1;
   21   22   0.7089   0.9373   0   250   250   250   0   0   1;
    3   23   0.4512   0.3083   0   250   250   250   0   0   1;
   23   24   0.8980   0.7091   0   250   250   250   0   0   1;
   24   25   0.8960   0.7011   0   250   250   250   0   0   1;
    6   26   0.2030   0.1034   0   250   250   250   0   0   1;
   26   27   0.2842   0.1447   0   250   250   250   0   0   1;
   27   28   0.6590   0.9337   0   250   250   250   0   0   1;
   28   29   0.8042   0.7006   0   250   250   250   0   0   1;
   29   30   0.5075   0.2585   0   250   250   250   0   0   1;
   30   31   0.9744   0.9630   0   250   250   250   0   0   1;
   31   32   0.3105   0.3619   0   250   250   250   0   0   1;
   32   33   0.3410   0.5302   0   250   250   250   0   0   1;
];
% Uncommenting the following line will create congestion in most situations
% mpc.branch(:, 6) = 1.6*ones(32, 1);
% Uncommenting the following line will eliminate congestion in most situations
mpc.branch(:, 6) = 250*ones(32, 1);
%% generator data
%	bus Pl   Pm    Ql    Qm      Ca     Cb
%   Objective function unit: $/(MW-h)
% mpc.gen = [
% 	1	0	0.9  -0.8  	 0.8	5.95    19.6;
%     2   0   1.1  -0.8  	 0.8	6.10    20.0;
%     3   0   1.1  -0.6  	 0.6	5.80    17.7;
%     6   0   0.8  -0.6  	 0.6	6.35    22.6;
%     9   0   0.8  -0.4  	 0.4	6.60    24.1;
% ];
mpc.gen = [
	1	0	20  -8  	 8	20;
    25  0 	0.66	0     0    0;
    28  0  0.6	0     0    0;
    31  0  0.67	0     0    0;
];
%% wind farm data
%	bus  cap	wf     sigma	cw
% mpc.wind = [
% 	2	 3.2    1.60    0.18    12.0;
%     3    2.2    1.10    0.11    12.0;
% ];