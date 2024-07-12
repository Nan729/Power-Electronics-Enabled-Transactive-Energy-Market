function BPF_SOCP = OPF_SOCP_new(Para)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_base = 10;			% kV
N_Bus = size(Para.bus,1);
% N_Gen = size(Para.gen,1);
% N_Wind= size(Para.wind,1);
N_Line= size(Para.branch,1); 
Ul = Para.bus(:,13).^2 * V_base^2; Um = Para.bus(:,12).^2 * V_base^2;	% pu -> kV^2
Pd = Para.bus(:,3);     Qd = Para.bus(:,4); 					% MW / MVAr
Gen_Pl=zeros(N_Bus,1);  Gen_Pl(Para.gen(:,1)) = Para.gen(:,2);	% MW / MVAr
Gen_Pm=zeros(N_Bus,1);  Gen_Pm(Para.gen(:,1)) = Para.gen(:,3);
Gen_Ql=zeros(N_Bus,1);  Gen_Ql(Para.gen(:,1)) = Para.gen(:,4);
Gen_Qm=zeros(N_Bus,1);  Gen_Qm(Para.gen(:,1)) = Para.gen(:,5);
Gen_Ca=zeros(N_Bus,1);  Gen_Ca(Para.gen(:,1)) = Para.gen(:,6);	% f: $/hr, P: MW
%Gen_Cb=zeros(N_Bus,1);  Gen_Cb(Para.gen(:,1)) = Para.gen(:,7);
Lr = Para.branch(:,3);  Lx = Para.branch(:,4); % Ohm
L_Cap = Para.branch(:,6); % Apply RateA - long-term capacity constraints
%cw=zeros(N_Bus,1);      cw(Para.wind(:,1)) = Para.wind(:,5);    
% MVP
% % Wind_l=zeros(N_Bus,1);  Wind_u=zeros(N_Bus,1);  Wind_u(Para.wind(:,1)) = Para.wind(:,3);	%MW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variable Declaration %%%%%%%%%%%%%%%%%%%%%%%%%%%
Bus_U = sdpvar(N_Bus,1);		% kV
Gen_P = sdpvar(N_Bus,1);		% MW
Gen_Q = sdpvar(N_Bus,1);		% MVAr
%Win_P = sdpvar(N_Bus,1);		% MW
Line_P = sdpvar(N_Line,1);		% MW
Line_Q = sdpvar(N_Line,1);		% MVAr
Line_I = sdpvar(N_Line,1);		% kA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% power flow model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Son_Line = find(Para.branch(:,1) == 1); % Only one line comes out of V-delta Bus
% Balancing Bus 1 Constraint (1a) - Active power balancing
Cons = (Line_P(Son_Line) == Gen_P(1));

% Other Buses Constraint (1a) - Active power balancing
for i = 2:N_Bus
Lout = find(Para.branch(:,1) == i);
Lin  = find(Para.branch(:,2) == i);  % For each bus, only one line comes in (tree)
% BNode_ID = Para.branch(Lin,1);
ENode_ID = Para.bus(i);    % Col i, Ln 1
if ~isempty(Lout)  % Regular line
    Cons = Cons + (sum(Line_P(Lout)) == Line_P(Lin) - Lr(Lin)*Line_I(Lin) + Gen_P(ENode_ID) - Pd(ENode_ID));
else               % Terminal line
    Cons = Cons + (0 == Line_P(Lin) - Lr(Lin)*Line_I(Lin) + Gen_P(ENode_ID)  - Pd(ENode_ID));
end
end

% Balancing Bus 1 Constraints (3) - Power cone and (1b) - Reactive power balancing
Cons = Cons + [ cone([2*Line_P(Son_Line); 2*Line_Q(Son_Line); Line_I(Son_Line) - Bus_U(1)], Line_I(Son_Line) + Bus_U(1));
                Gen_Q(1) == Line_Q(Son_Line)];

% Other Buses' Constraints
% (3) - Power cone, (1b) - Reactive power balancing, and (1c) - Voltage drop
for i = 2:N_Bus
Lout = find(Para.branch(:,1) == i);
Lin  = find(Para.branch(:,2) == i);
BNode_ID = Para.branch(Lin,1);
ENode_ID = Para.bus(i);
if ~isempty(Lout)  % Regular line
    Cons = Cons + (sum(Line_Q(Lout)) == Line_Q(Lin) - Lx(Lin)*Line_I(Lin) + Gen_Q(ENode_ID) - Qd(ENode_ID));  % Reactive Power Balacing
else               % Terminal line
    Cons = Cons + (Line_Q(Lin) - Lx(Lin)*Line_I(Lin) + Gen_Q(ENode_ID) - Qd(ENode_ID) == 0);
end
    Cons = Cons + [cone([2*Line_P(Lin); 2*Line_Q(Lin); Line_I(Lin) - Bus_U(BNode_ID)], Line_I(Lin) + Bus_U(BNode_ID)); % SOC
    Bus_U(ENode_ID) == Bus_U(BNode_ID) - 2*(Lr(Lin)*Line_P(Lin)+Lx(Lin)*Line_Q(Lin)) + (Lr(Lin)^2+Lx(Lin)^2)*Line_I(Lin)]; %voltage drop
end

% Cons-BND
Cons = Cons + [
        Gen_Pl <= Gen_P; Gen_P <= Gen_Pm; 
        Gen_Ql <= Gen_Q; Gen_Q <= Gen_Qm; 
        Ul <= Bus_U; Bus_U <=Um;   
        Line_I>=0;
        ];
 % Constraints on line flow
 % Sub the circumference with 4*N_Cir segments
 N_Cir = 32;
 Phi = (1:2:(2*N_Cir-1))'*90/2/N_Cir*pi/180;
 Phi0 = 90/2/N_Cir*pi/180;
 for i = 1:N_Line
     Cons = Cons + [sin(Phi)*Line_P(i) + cos(Phi)*Line_Q(i) <= repmat(L_Cap(i)*cos(Phi0), N_Cir, 1);
                    sin(Phi)*Line_P(i) - cos(Phi)*Line_Q(i) <= repmat(L_Cap(i)*cos(Phi0), N_Cir, 1);
                   -sin(Phi)*Line_P(i) + cos(Phi)*Line_Q(i) <= repmat(L_Cap(i)*cos(Phi0), N_Cir, 1);
                   -sin(Phi)*Line_P(i) - cos(Phi)*Line_Q(i) <= repmat(L_Cap(i)*cos(Phi0), N_Cir, 1);];
 end
 

%Cost = sum(Gen_Ca.*Gen_P.^2 + Gen_Cb.*Gen_P);

Cost = sum(Gen_Ca.*Gen_P);
ops = sdpsettings('verbose', 1, 'solver', 'mosek', 'verbose', 0);
sol = solvesdp(Cons,Cost,ops);
disp('SOCP optimization terminated.');
BPF_SOCP.Sol_Time = [sol.yalmiptime, sol.solvertime];

BPF_SOCP.Bus_U = double(Bus_U);
BPF_SOCP.Bus_V = sqrt(double(Bus_U));
BPF_SOCP.Gen_P = double(Gen_P(Para.gen(:,1)));
BPF_SOCP.Gen_Q = double(Gen_Q(Para.gen(:,1)));
%BPF_SOCP.Win_P = double(Win_P(Para.wind(:,1)));
BPF_SOCP.Line_P = double(Line_P);
BPF_SOCP.Line_Q = double(Line_Q);
BPF_SOCP.Line_I = double(Line_I);
BPF_SOCP.gap = zeros(N_Line,1);

for i = 1:N_Line
BNode_ID = Para.bus(:,1) == Para.branch(i,1);
% ENode_ID = find(Para.bus(:,1) == Para.branch(i,2));
BPF_SOCP.gap(i) = double(Lr(i)*(Line_I(i)*Bus_U(BNode_ID)-(Line_P(i).^2 + Line_Q(i).^2)));
end

BPF_SOCP.Cost = double(Cost);
BPF_SOCP.LMP = zeros(N_Bus,1);
for i = 1:N_Bus
BPF_SOCP.LMP(i) = dual(Cons(i));
end