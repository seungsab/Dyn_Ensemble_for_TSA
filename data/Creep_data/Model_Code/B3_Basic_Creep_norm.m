function J = B3_Basic_Creep_norm(x)
%% %%%%%%%%%%%%%%%% B3 CREEP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
method='basic';
%%%%%%%%%%%%%%%%%%%%%%%해석에 필요한 입력값 계산%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=t+t0;

wc=W/c;
sa=sand/(agg+sand);
w=c+W+agg+sand;

% COMPUPTE Volume-to-Surface ratio (mm)
switch SPEC
    case 0 % Cylinderical  Cross-Section
        volume=c_dia^2/4*pi*c_height;         %실험체부피
        surface=c_dia*pi*c_height;   %실험체표면적(위-아래는 제외)
        vs=volume/surface;           %부피/표면적
        
    case 1 % Square Cross-Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        volume=r_width*r_length*r_height;       %실험체부피
        surface=(r_width+r_length)*2*r_height;  %실험체표면적(위-아래는 제외)
        vs=volume/surface;                      %부피/표면적
        
end

switch method
    case 'basic'
        RH=100/100; 
    case 'drying'
        RH=RH/100;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ac=(agg+sand)/c;

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

e_s_inf=alpha1*alpha2*(1.9*10^(-2)*W^2.1*f28^-0.28+270);  %micro-strain

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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %회귀분석 변수
% switch method
%     case 'basic'
%         p1=0.6;
%         p2=185.4;
%         p3=0.29;
%         p4=20.3;
%     case 'drying'
%         p1=0.6;
%         p2=185.4;
%         p3=0.29;
%         p4=20.3;
%         p5=7.57;
%         p6=8.5;
% end

switch method
    case 'basic'
        p1=x(1);
        p2=x(2);
        p3=x(3);
        p4=x(4);
    case 'drying'
        p1=x(1);
        p2=x(2);
        p3=x(3);
        p4=x(4);
        p5=x(5);
        p6=x(6);          %8.5
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


q1=p1*10^6/E28;
q2=p2*c^0.5*f28^(-0.9);
q3=p3*(wc)^4*q2;
q4=p4*(ac)^(-0.7);

t0_load=max(t0,ts);

switch method
    case 'drying'
        D=2*vs/10;  %vs(cm)
        kt=p6*ts^(-0.08)*f28^(-0.25);
        tau_sh=kt*(ks*D)^2;
        e_sh_inf=e_s_inf*Elas(607,E28)/Elas(ts+tau_sh,E28);
        q5=p5*10^5/f28*(abs(e_sh_inf))^(-0.6);
end


for i=1:length(t)
    Q(i)=Qff(t0)*(1+(Qff(t0)/ZZ(t(i),t0))^gammaa(t0))^(-1/gammaa(t0));
end

switch method
    case 'basic'
        m=0.5;  n=0.1;
        for i=1:length(t)
            C0(i)=q2*Q(i)+q3*log(1+(t(i)-t0)^n)+q4*log(t(i)/t0);
            J_ALL(i)=q1+C0(i);
        end
    case 'drying'
        m=0.5;  n=0.1;
        for i=1:length(t)
            C0(i)=q2*Q(i)+q3*log(1+(t(i)-t0)^n)+q4*log(t(i)/t0);
            H(i)=1-(1-RH)*SS(t(i),ts,tau_sh);
            H_load=1-(1-RH)*SS(t0_load,ts,tau_sh);
            Cd(i)=q5*(exp(-8*H(i))-exp(-8*H_load))^0.5;
            J_ALL(i)=q1+C0(i)+Cd(i);
        end
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

function gamma=gammaa(t_loading)
gamma=1.7*(t_loading)^0.12+8;
end

function Qf=Qff(t_loading)
Qf=(0.086*(t_loading)^(2/9)+1.21*(t_loading)^(4/9))^(-1);
end

function S=SS(t,t0,tau_sh)
S=tanh(sqrt((t-t0)/tau_sh));
end

function Z=ZZ(t, t_loading)
m=0.5; 
n=0.1;
Z=(t_loading)^(-m)*log(1+(t-t_loading)^n);
end