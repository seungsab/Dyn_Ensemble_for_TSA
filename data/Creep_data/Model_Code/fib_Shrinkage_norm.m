function J = fib_Shrinkage_norm(x)
%% %%%%%%%%%%%%%%%% fib SHRINKAGE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
fcm=f28;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%회귀분석 변수
% switch c_type
%     case 2
%         alpha_bs=800;   
%     case 1
%         alpha_bs=700;   
%     case 3
%         alpha_bs=600;   
% end
% 
% fb=0.2;
% p1=220;
% fd=0.035;
% gamma=0.5;

fb=x(1);       % Normal range: 0.2 #1
alpha_bs=x(2); % Normal range: 600~800 #2
fd=x(3);       % Normal range: 0.035 #3
gamma=x(4);    % Normal range: 0.5 #4
p1=x(5);       % Normal range: 220 # 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch c_type
    case 2
        alpha_ds1=3;  alpha_ds2=0.013;
    case 1
        alpha_ds1=4;  alpha_ds2=0.012;
    case 3
        alpha_ds1=6;  alpha_ds2=0.012;
end

for i=1:length(t)
%%%%%%%%%%%%%%%%%%%자기수축%%%%%%%%%%%%%%%%%%%%%%%%
e_cbs0_fcm=-alpha_bs*((0.1*fcm)/(6+0.1*fcm))^2.5;
beta_bs_ts=1-exp(-fb*sqrt(ts));
beta_bs(i)=1-exp(-fb*sqrt(t(i)));
e_cbs(i)=e_cbs0_fcm*(beta_bs(i)-beta_bs_ts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%건조수축%%%%%%%%%%%%%%%%%%%%%%%
e_cds0_fcm=(p1+110*alpha_ds1)*exp(-alpha_ds2*fcm);

beta_s1=(35/fcm)^0.1;
if beta_s1<1
    beta_s1=1;
end

if 99<=RH
    beta_RH=0.25*beta_s1;                      %99<RH
else
    beta_RH=-1.55*(1-(RH/100)^3)*beta_s1;      %40<=RH<99
end

beta_ds(i)=((t(i)-ts)/(fd*(h)^2+(t(i)-ts)))^gamma;

e_cds(i)=e_cds0_fcm*beta_RH*beta_ds(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SH_ALL(i)=e_cbs(i)+e_cds(i);

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