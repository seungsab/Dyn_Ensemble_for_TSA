function J = ACI_Basic_Creep_norm(x)
%% %%%%%%%%%%%%%%%% ACI CREEP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - CALIBRATION INPUT: TUNING INPUT OF EMPIRICAL PREDICTION MODEL   ����! ���� ���� ��Ź
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
%%%%%%%%%%%%%%%%%%%%%%%�ؼ��� �ʿ��� �Է°� ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=t+t0;

wc=W/c;
sa=sand/(agg+sand);
w=c+W+agg+sand;

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
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%�������� ������ ���� ���(���� ������ �����Ҽ��� ũ������ ������)
if(cure==0)
    r_cure=1.25*(t0)^-0.118;  %���߾��(t0�� 7�̻�)
else
    r_cure=1.13*(t0)^-0.094;  %������(t0�� 1�̻�)
end

%�ܱ������ ���� ���
r_hum=1.27-0.0067*RH;  %�ܱ���� 40�ۼ�Ʈ �̻�

%����ü ���� ���� ���
%����-ǥ������
r_thick=2/3*(1+1.13*exp(-0.0213*vs));

%������
r_slump=0.82+0.00264*s;

%�ܰ����
r_agg=0.88+0.0024*sa;

%���ⷮ
r_air=0.46+0.09*a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ȸ�ͺм� ����(ȸ�ͺм� ������ ��� 3���� ������ �ͽ��ؼ� ��!)
% h=10;          % x(1);  % Normal range:6-30 days
% psi=0.6;       % x(2);  % Normal range:0.4-0.8  
% phi_inf=2.35;  % x(3);  % Normal range: 1.30-4.15
% Elas=Elastic(t0,w,f28, cure); %20,000-50,000

h=x(1);       % Normal range:6-30 days
psi=x(2);     % Normal range:0.4-0.8  
phi_inf=x(3); % Normal range: 1.30-4.15
Elas=x(4);    %20,000-50,000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%����ũ������� ���
phi_inf=phi_inf*r_cure*r_hum*r_thick*r_slump*r_agg*r_air;

%ũ�����Լ� ���
for i=1:length(t)
    phi(i)=(t(i)-t0)^psi/(h+(t(i)-t0)^psi)*phi_inf;         %ũ�������
%     J_ALL(i)=(1+phi(i))/Elastic(t0,w,f28, cure)*10^6;  %ũ�����Լ�
    J_ALL(i)=(1+phi(i))/Elas*10^6;  %ũ�����Լ�
end

J=J_ALL(:);
t=t-t0;

% plot(t,J_ALL);
% hold on
% plot(t, exp_data);

%% COMPUTE CREEP FUNCTION ON MEASUREMENT TIME
% INTERPOLATE THE CREEP (LINEAR INTERPOLATION)
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

function Ec=Elastic(t0, w, f28, cure)
if cure==0
    fc=(t0)/(4+0.85*(t0))*f28;
else
    fc=(t0)/(1+0.95*(t0))*f28;
end
Ec=0.043*(sqrt(w^3*fc));

end