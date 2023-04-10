function J = KCI_Drying_Creep_norm(x)
%% %%%%%%%%%%%%%%%% KCI CREEP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% COMPUPTE Volume-to-Surface ratio (mm)
switch SPEC
    case 0 % Cylinderical  Cross-Section
        Ac=c_dia^2/4*pi;         
        u=c_dia*pi;
        h=2*Ac/u;       
        
    case 1 % Square Cross-Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        Ac=r_width*r_length;  
        u=(r_width+r_length)*2;
        h=2*Ac/u;           
end

switch method
    case 'basic'
        RH=100; %외기습도(%) 
    case 'drying'
        RH=60; %외기습도(%) 
end

%%%%%%%%%탄성계수 계산%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Eci, Eci28]=Elastic(t0,w,f28);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%회귀분석 변수(회귀분석 변수는 요놈 3개를 적당히 믹스해서 씀!)
% bh=1.5*(1+(0.012*RH)^18)*h+250;
% if bh>=1500
%     bh=1500;
% end
% p1=16.8;  
% gamma=0.3;

p1=x(1);         % Normal range: 16.8
gamma=x(2);      % Normal range:0.17-0.43 
bh=x(3);         % Normal range:325-1500 days
Eci=x(4);        % 20,000-50,000
Eci28=x(4);      % 20,000-50,000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%크리프계수%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_dc_RH=1+(1-0.01*RH)/0.10/(h)^(1/3);  %외기습도
b_dc_fcu=p1/(f28)^0.5;                %강도
b_dc_t0=1/(0.1+t0^0.2);                %하중재하시점

Phi_0=b_dc_RH*b_dc_fcu*b_dc_t0;

for i=1:length(t)
b_c(i)=((t(i)-t0)/(bh+(t(i)-t0)))^gamma;
Phi(i)=b_c(i)*Phi_0;
end




%%%%%%%%%%크리프함수%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t)
    J_ALL(i)=(1/Eci+Phi(i)/Eci28)*10^6;  %크리프함수
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

function [Eci, Eci28]=Elastic(t0,mc,f28)

bcc=exp(0.35*(1-(28/t0)^0.5));
be=bcc^0.5;

Eci28=0.077*mc^1.5*(f28)^(1/3)*1.18;

Eci=Eci28*be;

end