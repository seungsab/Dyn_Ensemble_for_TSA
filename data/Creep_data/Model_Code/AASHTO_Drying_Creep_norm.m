function J = AASHTO_Drying_Creep_norm(x)
%% %%%%%%%%%%%%%%%% AASHTO CREEP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
[t,c_type,c,W,agg,sand,s,a,cure,RH,t0,ts,fck,f28,SPEC,c_dia,c_height,sh_width,sh_length,sh_height] =...
    deal(Input.t,Input.c_type,Input.cc,Input.water,Input.agg,Input.sand,Input.s,Input.a,Input.cure,Input.RH,Input.t0,...
    Input.ts,Input.fck,Input.f28,Input.SPEC.creep,Input.c_dia,Input.c_height,Input.sh_width,Input.sh_length,Input.sh_height);

x=Input.COEFF_nor(:,1).*x'+Input.COEFF_nor(:,2);


%% COMPUTE CREEP
method='drying';
%%%%%%%%%%%%%%%%%%%%%%%해석에 필요한 입력값 계산%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=t+t0;

wc=W/c;
sa=sand/(agg+sand);
w=c+W+agg+sand;

switch SPEC
    case 0 % Cylinderical  Cross-Section
        volume=c_dia^2/4*pi*c_height;         %실험체부피
        surface=c_dia*pi*c_height;   %실험체표면적(위-아래는 제외)
        vs=volume/surface/25.4;           %부피/표면적(inch)
        
    case 1 % Square Cross-Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        volume=r_width*r_length*r_height;       %실험체부피
        surface=(r_width+r_length)*2*r_height;  %실험체표면적(위-아래는 제외)
        vs=volume/surface/25.4;                      %부피/표면적(inch)       
end

switch method
    case 'basic'
        RH=100; %외기습도(%) 
    case 'drying'
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%콘크리트 강도 계산
f28=f28*0.145; %28일강도(ksi)

%%%%%%%%%%%%%%%%%%%%%보정계수%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%실험체 형상에 대한 계수
ks=1.45-0.13*vs;
if ks<1
    ks=1;
end

%외기습도에 대한 계수
khc=1.56-0.008*RH;

%콘크리트 강도에 대한 계수
if cure==0
    fci=t0/(4.0+0.85*t0)*f28; %강도계산(ACI 사용)
else
    fci=t0/(1.0+0.95*t0)*f28; %강도계산(ACI 사용)
end

kf=5/(1+fci);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%회귀분석 변수
% rate=61; %x(2);                   % 61
% p1=1.9;  %x(1);                   % 1.9
% Elas=Elastic(t0,w,f28);            %20,000-50,000

rate=x(1);                   % 61
p1=x(2);                   % 1.9
Elas=x(3);                 %20,000-50,000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%최종크리프계수 계산
phi_inf=p1*ks*khc*kf*t0^(-0.118);

%시간 함수 계산
for i=1:length(t)
    ktd(i)=(t(i)-t0)/(rate-4*fci+(t(i)-t0));
    %ktd(i)=(t(i)-t0)/(61-4*fci+(t(i)-t0));
end

%크리프함수 계수
for i=1:length(t)
    phi(i)=phi_inf*ktd(i);    %크리프계수
%     J_ALL(i)=(1+phi(i))/Elastic(t0,w,f28)*10^6;  %크리프함수
    J_ALL(i)=(1+phi(i))/Elas*10^6;  %크리프함수
end

J=J_ALL(:);
t=t-t0;

% plot(t,J_ALL);
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

function Ec=Elastic(t0, w, f28)
fc=(t0)/(4+0.85*(t0))*f28/0.145;
Ec=0.043*(sqrt(w^3*fc));
end