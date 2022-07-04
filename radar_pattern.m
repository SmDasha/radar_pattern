%% ������ ��������� �������������� �������� ������� (���)
clearvars
close all
clc
%% �������� ������
carrier_frequency= input('������� �������: '); % ������� �������
length= input('����� �����: ')/carrier_frequency; % ����� �����
sigma=input('��������� (��� ������): '); % ��������� (��� ������)
sigma_n=input('��������� (��� ����): '); % ��������� (��� ����)
num_elem_az = input('���������� ����������� � ������������ ���������: '); % ���������� ����������� � ������������ ���������
num_elem_el = input('���������� ����������� � ����������� ���������: '); % ���������� ����������� � ����������� ���������
num_elem = num_elem_az * num_elem_el;
k=input('����� ��������: '); % ����� ��������
elem_dist_az = input('���������� ����� ������������ � ������������ ��. � ������ ����: '); % ���������� ����� ������������ � ������������ ��.,
% �������� � ������ ����
elem_dist_el = input('���������� ����� ������������ � ����������� ��. � ������ ����'); % ���������� ����� ������������ � ����������� ��.,
% �������� � ������ ����
antenna_norm_dir = [input('����������� ������� � ������� �������, az: ');...
    input('����������� ������� � ������� �������, el: ')]; % ����������� ������� � ������� �������
% ������������ �������� ����������� �������������
k_u = db2pow(input('����������� �������� ������ ��������: ')); % ����������� �������� ������ ��������
az_signal=input('���� ������� ��������� �������: ')/180*pi; % ���� ������� ��������� �������
el_signal=input('���� ����� ��������� �������: ')/180*pi; % ���� ����� ��������� �������
scan_direction = [input('����������� ������������, az: '),...
    input('����������� ������������, el: ')]; % ����������� ������������ [���� �������;���� �����]
JAM_DIR = pi/180*[input('����������� ������, az: '), input('����������� ������, el: ');...
input('����������� ������, az: '), input('����������� ������, el: ')];% ����������� ������
%% ������� ��������� �������� �������
% ��������� ��������� ��������� � ������������ ��� ������������� ��. ������
% ������ � ������ �������.
elem_pos = zeros(2,num_elem_az*num_elem_el); % ��������� �������
elem_pos(1,:) = repmat(-(num_elem_az-1)*elem_dist_az/2:elem_dist_az:(num_elem_az-1)*elem_dist_az/2,...
1,num_elem_el); % ��������� ��������� � ������������ ���������
elem_pos(2,:) = reshape(repmat((num_elem_el-1)*elem_dist_el/2:-elem_dist_el:-(num_elem_el-1)*elem_dist_el/2,...
num_elem_az,1),1,[]); % ��������� ��������� � ������������ ���������
%% ������ ������� �������������
% ������� ���������
w_el = kaiser(num_elem_el,1);
w_az = kaiser(num_elem_az,1);
% w_el = hann(num_elem_el);
% w_az = hann(num_elem_az);
idx = 1;
k_scale = 2;
for i = 1:num_elem_az
    for j = 1:num_elem_el
        w(idx) = w_az(i)*w_el(j)/k_scale+(k_scale-1)/k_scale;
        idx = idx +1;
    end
end
w_hamming = w.';
% ������������ ��������� �������
V_signal=(elem_phase_scan_func([az_signal,el_signal],scan_direction,elem_pos)).';
V_signal=V_signal.*w_hamming;
% ������������ ������
V_p=(elem_phase_scan_func([JAM_DIR(:,1),JAM_DIR(:,2)],scan_direction,elem_pos)).';
%V_p=V_p.*w_hamming;
S_p=normrnd(0,sigma,numel(JAM_DIR)/2,k)+1i*normrnd(0,sigma,numel(JAM_DIR)/2,k);
SV_p=V_p*S_p;
% ������������ ����
randn('state',100);
z = normrnd(0,sigma_n,num_elem,k)+1i*normrnd(0,sigma_n,num_elem,k);
% ���+������
Y=SV_p+z;
% ���������� �������������� �������
R1=1/2*Y*Y';
R=R1+eye(num_elem)*5;
% ���������� �������� ������������
W=R\V_signal;
weights=sqrt(num_elem)*W/norm(W);
%%
% ������ �������� �� ��� ��� � ����������� �� ������ � �� ����������� ������������
u_weights=ones(num_elem,1);
gain_0 = abs(pattern_calc_func(scan_direction,u_weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
for idx = 1:numel(JAM_DIR)/2
    gain_j(idx) = abs(pattern_calc_func(JAM_DIR(idx,:),u_weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
end
disp('�� ���')
disp(['Antenna gain in scan direction: ', num2str(gain_0)])
disp(['Antenna gain in jammer direction: ', num2str(gain_j)])
% ������ �������� �� � ����������� �� ������ � �� ����������� ������������
gain_0 = abs(pattern_calc_func(scan_direction,weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
for idx = 1:numel(JAM_DIR)/2
    gain_j(idx) = abs(pattern_calc_func(JAM_DIR(idx,:),weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
end
disp('����� ���')
disp(['Antenna gain in scan direction: ', num2str(gain_0)])
disp(['Antenna gain in jammer direction: ', num2str(gain_j)])
disp(['Jam suppression: ', num2str(pow2db(gain_0./gain_j)), ' dB'])
angles_az = linspace(-pi/4,pi/4,180); % ����, ��� ������� ������������ ������
angles_el = linspace(-pi/4,pi/4,180); % ����, ��� ������� ������������ ������
pattern = zeros(numel(angles_az),numel(angles_el));
pattern_filtered = zeros(numel(angles_az),numel(angles_el));
for i = 1:numel(angles_az)
    for j = 1:numel(angles_el)
        angles = [angles_az(i), angles_el(j)];
        pattern(i,j) = abs(pattern_calc_func(angles,u_weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
        pattern_filtered(i,j) = abs(pattern_calc_func(angles,weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
    end
end
%% �������
% ������ ���������� ������ � ����������� "�������" � ����������� ������
imagesc(angles_az*180/pi,angles_el*180/pi,pattern) % ������� ������
title("���")
ylabel("�������� ��������, ����.")
xlabel("����, ����")
figure
imagesc(angles_az*180/pi,angles_el*180/pi,pattern_filtered) % ��������������� ������
title("��� ����� ���")
ylabel("�������� ��������, ����.")
xlabel("����, ����")
%����� ������� �������������
figure
imagesc(reshape(abs(weights),num_elem_az,num_elem_el))
%����� ������� � ������������ ���������
figure;
plot(angles_az*180/pi,pow2db(pattern(:,90)));hold on % ������� ������
plot(angles_az*180/pi,pow2db(pattern_filtered(:,90))); % ������� ������ � ��
ylim([-10,50])
grid on
title("��� � ������������ ���������")
ylabel("�������� ��������, ����.")
xlabel("����, ����")
xline(az_signal/pi*180,'g') % ������� �� ����������� ������������
legend ('��� ��� ���','��� ����� ��� � �������� kaiser','����������� �� ������������','����������� �� ������','����������� �� ������')
for i=1:numel(JAM_DIR(:,1)) % ������� ����������� �� ������
xline(JAM_DIR(i,1)/pi*180,'r')
end
% ����� ������� � ����������� ���������
figure;
plot(angles_el*180/pi,pow2db(pattern(90,:)));hold on % ������� ������
plot(angles_el*180/pi,pow2db(pattern_filtered(90,:))); % ������� ������ � ��
ylim([-10,50])
grid on
title("��� � ����������� ���������")
ylabel("�������� ��������, ����.")
xlabel("����, ����")
xline(el_signal/pi*180,'g') % ������� �� ����������� ������������
legend ('��� ��� ���','��� ����� ��� � �������� kaiser','����������� �� ������������')
for i=1:numel(JAM_DIR(:,2)) % ������� ����������� �� ������
    xline(JAM_DIR(i,2)/pi*180,'r')
end
%% ������ �������� ������
function[elem_phase_scan]=elem_phase_scan_func(angles,scan_direction,elem_pos)
elem_phase_scan = exp((sin(abs(angles-scan_direction))*elem_pos-...
fix(sin(abs(angles-scan_direction))*elem_pos))*2*pi*(-1i));
end
%% ������ ���
% ������ ������� ���� ���� �� ��������� ��� ����������� ������������,
% ������� ���� ���� - ����������, ������������� � ����� �����
function[pattern_calc]=pattern_calc_func(angles,weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u)
num_elem = num_elem_az*num_elem_el;
elem_pattern = k_u*num_elem*prod(sinc(angles*pi/180));
elem_phase_scan = exp((sin(abs(angles-scan_direction))-...
fix(sin(abs(angles-scan_direction))))*elem_pos*2*pi*(1i));
pattern_calc = sum(1/num_elem.*...
elem_pattern*elem_phase_scan*weights);
end