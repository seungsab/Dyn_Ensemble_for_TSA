function J = JSCE_Shrinkage_norm(x)
%% %%%%%%%%%%%%%%%% JSCE SHRINKAGE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%회귀분석 변수
% f=50;
% rate=0.108;
% power=0.56;

rate=x(1);   %0.108
power=x(2);  %0.56
f=x(3);      %50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%최종자기수축량 계산(고강도만 사용)%%%%%%%%%%%%%%%
e_as_inf=-3070*exp(-7.2*(W/c));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%최종 수축량 계산%%%%%%

if f28<=70
    e_sh_inf=-(-f+78*(1-exp(RH/100))+38*log(W)-5*(log(volume/surface/10))^2)*10;
else
    nu=10^-4*(15*exp(0.007*f28)+0.25*W);
    if c_type==1
        a=11;   % ordinary, low heat =11, high-early-strength =15
    elseif c_type==3
        a=15;
    else
        a=11;
    end
    e_ds_ro=(a*(1-RH/100)*W)/(1+150*exp(-500/f28));
    e_ds_inf=-e_ds_ro/(1+nu*ts);
end


%%%%%%%%%%수축량 계산%%%%%%%%%%%
if f28<=70
    for i=1:length(t)
        SH_ALL(i)=(1-exp(-rate*(t(i)-ts)^power))*e_sh_inf;
    end
else
    if 0.2<=(W/c)<0.23
        aa=(1.5-1.2)/(0.23-0.2)*((W/c)-0.2)+1.2;
        bb=0.4;
    elseif 0.23<=(W/c)<0.3
        aa=(0.6-1.5)/(0.3-0.23)*((W/c)-0.23)+1.5;
        bb=(0.5-0.4)/(0.3-0.23)*((W/c)-0.23)+0.4;
    elseif 0.3<=(W/c)<0.4
        aa=(0.1-0.6)/(0.4-0.3)*((W/c)-0.3)+0.6;
        bb=(0.7-0.5)/(0.4-0.3)*((W/c)-0.3)+0.5;
    elseif 0.4<=(W/c)<0.5
        aa=(0.03-0.1)/(0.5-0.4)*((W/c)-0.4)+0.1;
        bb=(0.8-0.7)/(0.5-0.4)*((W/c)-0.4)+0.7;
    else
        aa=0.03;
        bb=0.8;
    end
    
     
    beta=(4*W*sqrt(vs))/(100+0.7*ts);
    for i=1:length(t)
        e_sh(i)=(e_ds_inf*(t(i)-ts))/(beta+(t(i)-ts));
        e_ast(i)=e_as_inf*((1-exp(-aa*(t(i))^bb))-(1-exp(-aa*(ts)^bb)));
        SH_ALL(i)=e_sh(i)+e_ast(i);
    end
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