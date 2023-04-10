function J = ACI_Shrinkage_norm(x)
%% %%%%%%%%%%%%%%%% ACI CREEP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        vs=volume/surface;           %부피/표면적(inch)
    case 1
        volume=sh_width*sh_length*sh_height;       %실험체부피
        surface=(sh_width+sh_height)*2*sh_length;  %실험체표면적(위-아래는 제외)
        vs=volume/surface;                      %부피/표면적(inch)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%외기습도에 대한 계수(%)
if RH<=80
    r_hum=1.40-0.0102*RH;  %40<=RH<=80
else
    r_hum=3.00-0.030*RH;   %80<RH<=100
end

%실험체 형상에 대한 계수(mm)
%부피-표면적비
r_thick=1.2*exp(-0.00472*(vs));


%슬럼프(mm)
r_slump=0.89+0.00161*s;

%잔골재비(%)
if sa<=50
    r_agg=0.30+0.014*sa;    %sa<=50
else
    r_agg=0.90+0.002*sa;    %sa>50
end

%시멘트양(kg/m3)
r_c=0.75+0.00061*c;

%공기량(%)
r_air=0.95+0.008*a;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%회귀분석 변수
% if(cure==0)
%     f=35;  %수중양생(t0는 7이상)
% else
%     f=55;  %증기양생(t0는 1이상)
% end
% 
% a=1;   %x(2); %Normal range:0.9-1.10 
% shu0=780;  %x(3); %Normal range:415-1070

f=x(1);  %Normal range:20-130
a=x(2); %Normal range:0.9-1.10 
shu0=x(3); %Normal range:415-1070
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%최종 수축량 계산%%%%%%
e_shu=-shu0*r_hum*r_thick*r_slump*r_agg*r_c*r_air;


%%%%%%%%%%수축량 계산%%%%
for i=1:length(t)
    SH_ALL(i)=(t(i)-ts)^a/(f+(t(i)-ts)^a)*e_shu;
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
% if strcmp(Input.Infer,'GA')
%     J = repmat(J,Input.N_Y,1); J = sum((Input.Y - J).^2);
% elseif strcmp(Input.Infer,'DREAM')
%     J = repmat(J,Input.N_Y,1);
%end
end