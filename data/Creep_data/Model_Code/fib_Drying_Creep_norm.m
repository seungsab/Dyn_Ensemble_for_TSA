function J = fib_Drying_Creep_norm(x)
%% %%%%%%%%%%%%%%%% fib CREEP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - CALIBRATION INPUT: TUNING INPUT OF EMPIRICAL PREDICTION MODEL   ����! ���� ���� ��Ź
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
%%%%%%%%%%%%%%%%%%%%%%%�ؼ��� �ʿ��� �Է°� ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=t+t0;

wc=W/c;
sa=sand/(agg+sand);
w=c+W+agg+sand;
fcm=f28;

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
        RH=100; %�ܱ����(%) 
    case 'drying'
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_fcm=(35/fcm)^0.5;

%ź����� ���
[Eci, Eci28]=Elastic(t0,fcm); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ȸ�ͺм� ����(ȸ�ͺм� ������ ��� 3���� ������ �ͽ��ؼ� ��!)
bh=1.5*h+250*alpha_fcm;
if bh>=1500*alpha_fcm
    bh=1500*alpha_fcm;
end
gamma=1/(2.3+(3.5/sqrt(t0)));

% p1=30;
% p2=1.8;
% p3=412;

switch method
    case 'basic'
        p2=x(1);       % 1.8 #1
        p1=x(2);       % 30 #2
        p3=412;       % 412 #3       
        Eci=x(3);      % 20,000-50,000
        Eci28=x(3);    % 20,000-50,000
    case 'drying'
        p2=x(1);       % 1.8 #1
        p1=x(2);       % 30 #2
        p3=x(3);       % 412 #3
        bh=x(4);       % Normal range:200-2000 days #4
        gamma=x(5);    % Normal range:0.17-0.43 #5
        Eci=x(6);      % 20,000-50,000
        Eci28=x(6);    % 20,000-50,000
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%�⺻ũ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_bc_fcm=p2/(fcm)^0.7;

%%%%%%%%%%%%����ũ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_dc_fcm=p3/(fcm)^1.4;                %����
b_dc_RH=(1-RH/100)/(0.1*h/100)^(1/3);  %�ܱ����
b_dc_t0=1/(0.1+t0^0.2);                %�������Ͻ���

for i=1:length(t)
    %%%%%%%%%%%%%ũ���� �ð��Լ� ���%%%%%%%%%%%%%%%%%%%%%
    b_bc_t(i)=log((p1/t0+0.035)^2*(t(i)-t0)+1);  %�⺻
    b_dc_t(i)=((t(i)-t0)/(bh+(t(i)-t0)))^gamma;     %����
    
    %%%%%%%%%%%%%%ũ�������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Phi_bc(i)=b_bc_fcm*b_bc_t(i);  %�⺻ũ����
    Phi_dc(i)=b_dc_fcm*b_dc_RH*b_dc_t0*b_dc_t(i);  %����ũ����
    
    %%%���� ũ��������� �̿� ���������� ũ������� ���%%%%%
    switch method
        case 'basic'
            Phi(i)=Phi_bc(i);
        case 'drying'
            Phi(i)=Phi_bc(i)+Phi_dc(i);
    end
end


%ũ�����Լ�, ���
for i=1:length(t)
    J_ALL(i)=(1/Eci+Phi(i)/Eci28)*10^6;  %ũ�����Լ�
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

function [Eci, Eci28]=Elastic(t0,f28)

bcc=exp(0.25*(1-(28/t0)^0.5));
be=bcc^0.5;

alphaE=1.0; %���� ������ ���� ��, 0.7-1.2
Eci28=21.5*10^3*alphaE*(f28/10)^(1/3);

Eci=Eci28*be;

end