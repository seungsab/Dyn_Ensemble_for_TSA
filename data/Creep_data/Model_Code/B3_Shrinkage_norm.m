function J = B3_Shrinkage_norm(x)
%% %%%%%%%%%%%%%%%% B3 SHRINKAGE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - CALIBRATION INPUT: TUNING INPUT OF EMPIRICAL PREDICTION MODEL   차씨! 변수 정의 부탁
% h            :
% psi          :
% psi_inf      :
%% - VARIABLE INPUT
% Input.method : BASIC or DRYING CREEP
% Input.t      : time (d)
% Input.Y      : Experimental data [Matrix = (#point) X (n_rep+1)], First column is the time

%Input.c_type  : Cement classification (1, 2, 3, 5)
%Input.cc      : Cement content (kg/m3)
%Input.water   : Water content (kg/m3)
%Input.agg     : Aggregate content (kg/m3)
%%%Input.sand    : Sand content (kg/m3)
%Input.s       : Slump (mm)
%Input.a       : Air content (%)

% Input.cure   : Curing condition % 0 is mositure curing, except 0 is steam curing
% Input.RH     : Relative humidity (%)
% Input.t0     : Loading age (d)
% Input.ts     : Exposure time to air (day)

%Input.fck     : Specific compressive strength (MPa)
% Input.f28    : Compressive strength of standard cylinders (at 28 days)

% Input.SPEC   : Concrete Specimens (0: Cylinder // 1: Square Cross-Section)
%Input.c_dia   : Cylinder diameter  (mm)
%Input.c_height: Cylinder height    (mm)
%Input.r_width : Rectangle(square) width     (mm)
%Input.r_length: Rectangle(square) length    (mm)
%Input.r_height: Rectangle(square) height    (mm)

%Input.shc_dia   : Cylinder diameter  (mm)
%Input.shc_height: Cylinder height    (mm)
%Input.sh_width : Rectangle(square) width     (mm)
%Input.sh_length: Rectangle(square) length    (mm)
%Input.sh_height: Rectangle(square) height    (mm)
%% - OUTPUT
% J            : Creep function 
% Phi          : Creep coefficient
%% HISTORY
% CODED BY SS (20170308) BASED ON Cha's CODE
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Input

%% ASSIGN VARIALBE INPUTS
Input.method = 'a';
[t,c_type,c,W,agg,sand,s,a,cure,RH,t0,ts,fck,f28,SPEC,shc_dia,shc_height,sh_width,sh_length,sh_height] =...
    deal(Input.t,Input.c_type,Input.cc,Input.water,Input.agg,Input.sand,Input.s,Input.a,Input.cure,Input.RH,Input.t0,...
    Input.ts,Input.fck,Input.f28,Input.SPEC.shrink,Input.c_dia,Input.c_height,Input.sh_width,Input.sh_length,Input.sh_height);

x=Input.COEFF_nor(:,1).*x'+Input.COEFF_nor(:,2);

%% COMPUTE SHRINKAGE
%%%%%%%%%%%%%%%%%%%%%%%해석에 필요한 입력값 계산%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=t+ts;

wc=W/c;
sa=sand/(agg+sand);
w=c+W+agg+sand;

% COMPUPTE Volume-to-Surface ratio (mm)
switch SPEC
    case 0
        volume=shc_dia^2/4*pi*shc_height;         %실험체부피
        surface=shc_dia*pi*shc_height;   %실험체표면적(위-아래는 제외)
        vs=volume/surface;           %부피/표면적
    case 1
        volume=sh_width*sh_length*sh_height;       %실험체부피
        surface=(sh_width+sh_height)*2*sh_length;  %실험체표면적(위-아래는 제외)
        vs=volume/surface;                      %부피/표면적
end

RH=RH/100;

ac=(agg+sand)/c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Basic creep
E28=4734*f28^0.5;  


%Drying creep
switch c_type
    case 1
        alpha1=1;
    case 2
        alpha1=0.85;
    case 3
        alpha1=1.1;
end

switch cure
    case 1
        alpha2=0.75;
    case '밀봉또는기중'
        alpha2=1.2;
    case 0
        alpha2=1.0;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %회귀분석 변수
% power=0.5;
% p2=8.5;
% p1=270;


power=x(1);
p2=x(2);
p1=x(3);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


e_s_inf=alpha1*alpha2*(1.9*10^(-2)*W^2.1*f28^-0.28+p1);  %micro-strain

if SPEC==0
    specimen='infinite_cylinder';
else
    specimen='infinite_square_prism';
end


switch specimen
    case 'infinite_slab'
        ks=1;
    case 'infinite_cylinder'
        ks=1.15;
    case 'infinite_square_prism'
        ks=1.25;
    case 'sphere'
        ks=1.3;
    case 'cube'
        ks=1.55;
end

if RH<=0.98
    k_h=1-RH^3;
else
    k_h=(-0.2-0.0588)/(1-0.98)*(RH-1)-0.2;
end

t0_load=max(t0,ts);
D=2*vs/10;  %vs(cm)
kt=p2*ts^(-0.08)*f28^(-0.25);

tau_sh=kt*(ks*D)^2;

e_sh_inf=-e_s_inf*Elas(607,E28)/Elas(ts+tau_sh,E28);

for i=1:length(t)
    SH_ALL(i)=e_sh_inf*k_h*tanh(((t(i)-ts)/tau_sh)^power);
end

J=SH_ALL(:);
t=t-ts;

% plot(t,SH_ALL);
% hold on
% plot(t, exp_data);

% %% COMPUTE CREEP FUNCTION ON MEASUREMENT TIME
% % INTERPOLATE THE CREEP (LINEAR INTERPOLATION)
% J = interp1(t,SH_ALL,Input.t);
% if any(isnan(J))
%     J = SH_ALL';
% end
% 
% 
% if strcmp(Input.Infer,'GA')
%     J = repmat(J,Input.N_Y,1); J = sum((Input.Y - J).^2);
% elseif strcmp(Input.Infer,'DREAM')
%     J = repmat(J,Input.N_Y,1);
% end

end

function E=Elas(t,E28)
E=E28*(t/(4+0.85*t))^0.5;
end