function J = AASHTO_Shrinkage_norm(x)
%% %%%%%%%%%%%%%%%% AASHTO SHRINKAGE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% COMPUPTE Volume-to-Surface ratio (mm)
t=t+ts;

wc=W/c;
sa=sand/(agg+sand);
w=c+W+agg+sand;

switch SPEC
    case 0
        volume=shc_dia^2/4*pi*shc_height;         %실험체부피
        surface=shc_dia*pi*shc_height;   %실험체표면적(위-아래는 제외)
        vs=volume/surface/25.4;           %부피/표면적(inch)
    case 1
        volume=sh_width*sh_length*sh_height;       %실험체부피
        surface=(sh_width+sh_height)*2*sh_length;  %실험체표면적(위-아래는 제외)
        vs=volume/surface/25.4;                      %부피/표면적(inch)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%콘크리트 강도
f28=f28*0.145; %28일강도(ksi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%보정계수%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%외기습도에 대한 계수
k_hs=2.00-0.014*RH;

%콘크리트 강도에 대한 계수
if cure==0
    fci=ts/(4.0+0.85*ts)*f28; %강도계산(ACI 사용)
else
    fci=ts/(1.0+0.95*ts)*f28; %강도계산(ACI 사용)
end

k_f=5/(1+fci);

%실험체 형상에 대한 계수
k_s=1.45-0.13*vs;
if k_s<1
    k_s=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%회귀분석 변수
% esh0=0.48; %x(1);  % 0.48
% rate=61;   %x(2);  % 61

rate=x(1);  % 61
esh0=x(2);  % 0.48
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%최종수축량 계산
e_sh_inf=-esh0*10^-3*k_hs*k_f*k_s*10^6;

%시간 함수 및 수축 계산
for i=1:length(t)
    %k_s(i)=(t(i)/(26*exp(0.36*vs)+t(i)))/(t(i)/(45+t(i)))*(1064-94*vs)/923;
    k_td(i)=(t(i)-ts)/(rate-4*fci+(t(i)-ts));
    SH_ALL(i)=e_sh_inf*k_td(i);
end

J=SH_ALL(:);
t=t-ts;

% plot(t,SH_ALL);
% hold on
% plot(t, exp_data);

% %% COMPUTE CREEP FUNCTION ON MEASUREMENT TIME
% % INTERPOLATE THE CREEP (LINEAR INTERPOLATION)
% J = interp1(t,J_ALL,Input.t);
% if any(isnan(J))
%     J = J_ALL';
% end
% 
% if strcmp(Input.Infer,'GA')
%     J = repmat(J,Input.N_Y,1); J = sum((Input.Y - J).^2);
% elseif strcmp(Input.Infer,'DREAM')
%     J = repmat(J,Input.N_Y,1);
% end
end