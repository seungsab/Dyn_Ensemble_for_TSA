function J = GL_Basic_Creep_norm(x)
%% %%%%%%%%%%%%%%%% GL CREEP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        RH=96/100; %외기습도(%) 
    case 'drying'
        RH=RH/100;
end


switch c_type
    case 1 
        ss=0.335; k=1.0;
    case 2 
        ss=0.4; k=0.75;
    case 3
        ss=0.13; k=1.15;
end

[E_cmt, E_cm28]=Elastic(t0,f28,ss); %탄성계수 계산
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%회귀분석 변수
% p1=2;
% p2=0.3;
% p3=14;
% p4=2.5;
% p5=0.12;
% p6=0.5;
% p7=E_cmt                 %20,000-50,000
% p8=E_cm28                %20,000-50,000

switch method
    case 'basic'
        p1=x(1);
        p2=x(2);
        p3=x(3);
        p4=2.5;
        p5=0.12;
        p6=0.5;
        p7=x(4);     %E_cmt   %20,000-50,000
        p8=x(4);     %E_cm28  %20,000-50,000
    case 'drying'
        p1=x(1);
        p2=x(2);
        p3=x(3);
        p4=x(4);
        p5=x(5);
        p6=x(6);
        p7=x(7);     %E_cmt   %20,000-50,000
        p8=x(7);     %E_cm28  %20,000-50,000
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%크리프계수%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_tc=(1-((t0-ts)/((t0-ts)+p5*(volume/surface)^2))^0.5)^0.5;

for i=1:length(t)
    Phi_28(i)=Phi_tc*(p1*(t(i)-t0)^p2/((t(i)-t0)^p2+p3)+(7/t0)^0.5*((t(i)-t0)/((t(i)-t0)+7))^0.5...
        +p4*(1-1.086*RH^2)*((t(i)-t0)/((t(i)-t0)+p5*(volume/surface)^2))^p6);
end


%%%%%%%%%%%크리프함수%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(t)
    J_ALL(i)=(1/E_cmt+Phi_28(i)/E_cm28)*10^6;  %크리프함수
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

function [E_cmt,E_cm28]=Elastic(t0,f28,ss)

%f_cm28=1.1*fc+5.0;   %fc: 설계기준압축강도(MPa)
f_cm28=f28;

beta_e=exp(ss/2*(1-sqrt(28/t0)));

f_cmt=(beta_e)^2*f_cm28;   

E_cmt=3500+4300*sqrt(f_cmt); %MPa
E_cm28=3500+4300*sqrt(f_cm28); %MPa

end