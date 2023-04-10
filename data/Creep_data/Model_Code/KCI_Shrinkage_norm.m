function J = KCI_Shrinkage_norm(x)
%% %%%%%%%%%%%%%%%% KCI SHRINKAGE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - CALIBRATION INPUT: TUNING INPUT OF EMPIRICAL PREDICTION MODEL   차씨! 변수 정의 부탁
% bh            :
% gamma         :
% Phi           :
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
fcu=f28;

% COMPUPTE Volume-to-Surface ratio (mm)
switch SPEC
    case 0
        Ac=shc_dia^2*pi/4;
        u=shc_dia*pi;
        h=2*Ac/u;       
    case 1       
        Ac=sh_width*sh_height;
        u=(sh_width+sh_height)*2;
        h=2*Ac/u;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%회귀분석 변수(회귀분석 변수는 요놈 3개를 적당히 믹스해서 씀!)
% p1=160;
% power=0.5;
% f=0.035;

f=x(1);          % Normal range:   0.035
power=x(2);      % Normal range:    0.5
p1=x(3);         % Normal range:   160
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch c_type
    case 1 || 5
        beta_sc=5;
    case 2
        beta_sc=4;
    case 3
        beta_sc=8;
end

e_s_fcu=p1+10*beta_sc*(9-fcu/10);

if RH<99
    beta_RH=-1.55*(1-(RH/100)^3);
else
    beta_RH=0.25;
end

%%%%%%%%%%%%%%%%%%최종수축량계산%%%%%%%%%%%%%%%%%%%%%%%%%
e_sho=e_s_fcu*beta_RH;

%%%%%%%%%%%%%%%%%%%%%%%%%수축량계산%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t)
    beta_s(i)=((t(i)-ts)/(f*(h)^2+(t(i)-ts)))^power;
    SH_ALL(i)=e_sho*beta_s(i);
end

J=SH_ALL(:);
t=t-ts;

% plot(t,SH_ALL);
% hold on
% plot(t, exp_data);

% %% COMPUTE SHRINKAGE ON MEASUREMENT TIME
% % INTERPOLATE THE CREEP (LINEAR INTERPOLATION)
% J = interp1(t,SH_ALL,Input.t);
% if any(isnan(J))
%     J = SH_ALL';
% end
% 
% if strcmp(Input.Infer,'GA')
%     J = repmat(J,Input.N_Y,1); J = sum((Input.Y - J).^2);
% elseif strcmp(Input.Infer,'DREAM')
%     J = repmat(J,Input.N_Y,1);
% end
end