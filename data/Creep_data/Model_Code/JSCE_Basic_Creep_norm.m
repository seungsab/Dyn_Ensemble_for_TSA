function J = JSCE_Basic_Creep_norm(x)
%% %%%%%%%%%%%%%%%% JSCE CREEP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
method='basic';
%%%%%%%%%%%%%%%%%%%%%%%�ؼ��� �ʿ��� �Է°� ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=t+t0;

wc=W/c;
sa=sand/(agg+sand);
w=c+W+agg+sand;

% COMPUPTE Volume-to-Surface ratio (mm)
switch SPEC
    case 0 % Cylinderical  Cross-Section
        volume=c_dia^2/4*pi*c_height;         %����ü����
        surface=c_dia*pi*c_height;   %����üǥ����(��-�Ʒ��� ����)
        vs=volume/surface;           %����/ǥ����
        
    case 1 % Square Cross-Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        volume=r_width*r_length*r_height;       %����ü����
        surface=(r_width+r_length)*2*r_height;  %����üǥ����(��-�Ʒ��� ����)
        vs=volume/surface;                      %����/ǥ����       
end

switch method
    case 'basic'
        RH=100; %�ܱ����(%) 
    case 'drying'
                %�ܱ����(%) 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ũ��Ʈ ������ ���� ���
if cure==0
    fc_loading=t0/(4.0+0.85*t0)*f28; %�������(ACI ���)
else
    fc_loading=t0/(1.0+0.95*t0)*f28; %�������(ACI ���)
end

%���� ������ ��ũ��Ʈ ����   55MPa����(70���� ��밡��), ������ ���� �Լ��� ���� ���� 55�̻� 80MPa ����
%fc_loading=55;
%�����ø�Ʈ��(kg/m3) 260<=c<=500
%��-�ø�Ʈ ��  ������ 0.4<=w/c<=0.65, ������ ���� ����

%% COMPUTE CREEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ȸ�ͺм� ����(���������� ���� 55MPa(70) ����)
% switch method
%     case 'basic'
%         f_b=15;
%         rate=0.09;
%         power=0.6;
%         [Ects]=Elastic(fc_loading); 
%
%         phi_bc_inf=(f_b*(c+W)^2.0*(wc)^2.4*(log(t0))^-0.67)*10^-10*10^6;
%        
%         
%     case 'drying'
%         f_b=15;
%         rate=0.09;
%         power=0.6;
%         f_d=4500;
%         [Ects]=Elastic(fc_loading);

%         phi_bc_inf=(f_b*(c+W)^2.0*(wc)^2.4*(log(t0))^-0.67)*10^-10*10^6;
%         phi_dc_inf=(f_d*(c+W)^1.4*(wc)^4.2*(log(volume/surface/10))^-2.2*(1-RH/100)^0.36*(ts)^-0.30)*10^-10*10^6;
% end


switch method
    case 'basic'
        rate=x(1);   %Normal range: 0.09
        power=x(2);   %% 0.6
        f_b=x(3);    %Normal range: 15
        Ects=x(4);   %20,000-50,000
        
        phi_bc_inf=(f_b*(c+W)^2.0*(wc)^2.4*(log(t0))^-0.67)*10^-10*10^6;
        
    case 'drying'
        rate=x(1);   %Normal range: 0.09
        power=x(2);   %% 0.6
        f_b=x(3);    %Normal range: 15                
        f_d=x(4);    %Normal range: 4500
        Ects=x(5);   %20,000-50,000
        
        phi_bc_inf=(f_b*(c+W)^2.0*(wc)^2.4*(log(t0))^-0.67)*10^-10*10^6;
        phi_dc_inf=(f_d*(c+W)^1.4*(wc)^4.2*(log(volume/surface/10))^-2.2*(1-RH/100)^0.36*(ts)^-0.30)*10^-10*10^6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%����ũ������� ���%%%%%%%%%%%%%%

switch method
    case 'basic'
        phi_inf=phi_bc_inf;
    case 'drying'
        phi_inf=phi_bc_inf+phi_dc_inf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ũ�����Լ� ���
if fc_loading<=70
    for i=1:length(t)
        J_ALL(i)=(1-exp(-rate*(t(i)-t0)^power))*phi_inf+1/Ects*10^6;
    end 
else
    for i=1:length(t)
        J_ALL(i)=(4*W*(1-RH/100)+350)/(12+fc_loading)*log(t(i)-t0+1)+1/Ects*10^6;
    end
end

J=J_ALL(:);
t=t-t0;
% 
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

function [Ects]=Elastic(fc_loading)

Ects=(41.6-39.2*0.973^fc_loading)*1000;

end