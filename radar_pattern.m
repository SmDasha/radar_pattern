%% Расчёт диаграммы направленности антенной решётки (ДНА)
clearvars
close all
clc
%% Исходные данные
carrier_frequency= input('несущая частота: '); % несущая частота
length= input('длина волны: ')/carrier_frequency; % длина волны
sigma=input('дисперсия (для помехи): '); % дисперсия (для помехи)
sigma_n=input('дисперсия (для шума): '); % дисперсия (для шума)
num_elem_az = input('количество излучателей в азимутальной плоскости: '); % Количество излучателей в азимутальной плоскости
num_elem_el = input('количество излучателей в угломестной плоскости: '); % Количество излучателей в угломестной плоскости
num_elem = num_elem_az * num_elem_el;
k=input('число отсчетов: '); % число отсчетов
elem_dist_az = input('расстояние между излучателями в азимутальной пл. в длинах волн: '); % Расстояние между излучателями в азимутальной пл.,
% заданное в длинах волн
elem_dist_el = input('расстояние между излучателями в угломестной пл. в длинах волн'); % Расстояние между излучателями в угломестной пл.,
% заданное в длинах волн
antenna_norm_dir = [input('направление нормали к полотну антенны, az: ');...
    input('направление нормали к полотну антенны, el: ')]; % Направление нормали к полотну антенны
% относительно нулевого направления радиолокатора
k_u = db2pow(input('коэффициент усиления одного элемента: ')); % коэффициент усиления одного элемента
az_signal=input('угол азимута полезного сигнала: ')/180*pi; % угол азимута полезного сигнала
el_signal=input('угол места полезного сигнала: ')/180*pi; % угол места полезного сигнала
scan_direction = [input('направление сканирования, az: '),...
    input('направление сканирования, el: ')]; % направление сканирования [угол азимута;угол места]
JAM_DIR = pi/180*[input('направление помехи, az: '), input('направление помехи, el: ');...
input('направление помехи, az: '), input('направление помехи, el: ')];% направление помехи
%% Задание геометрии антенной решётки
% Вычисляем положение элементов в пространстве для прямоугольной АР. Отсчёт
% ведётся с центра решётки.
elem_pos = zeros(2,num_elem_az*num_elem_el); % аллокация массива
elem_pos(1,:) = repmat(-(num_elem_az-1)*elem_dist_az/2:elem_dist_az:(num_elem_az-1)*elem_dist_az/2,...
1,num_elem_el); % положения элементов в азимутальной плоскости
elem_pos(2,:) = reshape(repmat((num_elem_el-1)*elem_dist_el/2:-elem_dist_el:-(num_elem_el-1)*elem_dist_el/2,...
num_elem_az,1),1,[]); % положения элементов в азимутальной плоскости
%% Расчет весовых коэффициентов
% Оконная обработка
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
% Формирование полезного сигнала
V_signal=(elem_phase_scan_func([az_signal,el_signal],scan_direction,elem_pos)).';
V_signal=V_signal.*w_hamming;
% Формирование помехи
V_p=(elem_phase_scan_func([JAM_DIR(:,1),JAM_DIR(:,2)],scan_direction,elem_pos)).';
%V_p=V_p.*w_hamming;
S_p=normrnd(0,sigma,numel(JAM_DIR)/2,k)+1i*normrnd(0,sigma,numel(JAM_DIR)/2,k);
SV_p=V_p*S_p;
% Формирование шума
randn('state',100);
z = normrnd(0,sigma_n,num_elem,k)+1i*normrnd(0,sigma_n,num_elem,k);
% Шум+помеха
Y=SV_p+z;
% Нахождение корреляционной матрицы
R1=1/2*Y*Y';
R=R1+eye(num_elem)*5;
% Нахождение весового коэффициента
W=R\V_signal;
weights=sqrt(num_elem)*W/norm(W);
%%
% Оценка значений ДН без АПФ в направлении на помехи и по направлению сканирования
u_weights=ones(num_elem,1);
gain_0 = abs(pattern_calc_func(scan_direction,u_weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
for idx = 1:numel(JAM_DIR)/2
    gain_j(idx) = abs(pattern_calc_func(JAM_DIR(idx,:),u_weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
end
disp('До АПФ')
disp(['Antenna gain in scan direction: ', num2str(gain_0)])
disp(['Antenna gain in jammer direction: ', num2str(gain_j)])
% Оценка значений ДН в направлении на помехи и по направлению сканирования
gain_0 = abs(pattern_calc_func(scan_direction,weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
for idx = 1:numel(JAM_DIR)/2
    gain_j(idx) = abs(pattern_calc_func(JAM_DIR(idx,:),weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
end
disp('После АПФ')
disp(['Antenna gain in scan direction: ', num2str(gain_0)])
disp(['Antenna gain in jammer direction: ', num2str(gain_j)])
disp(['Jam suppression: ', num2str(pow2db(gain_0./gain_j)), ' dB'])
angles_az = linspace(-pi/4,pi/4,180); % Углы, для которых производится расчёт
angles_el = linspace(-pi/4,pi/4,180); % Углы, для которых производится расчёт
pattern = zeros(numel(angles_az),numel(angles_el));
pattern_filtered = zeros(numel(angles_az),numel(angles_el));
for i = 1:numel(angles_az)
    for j = 1:numel(angles_el)
        angles = [angles_az(i), angles_el(j)];
        pattern(i,j) = abs(pattern_calc_func(angles,u_weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
        pattern_filtered(i,j) = abs(pattern_calc_func(angles,weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u));
    end
end
%% Графики
% Задача подавления помехи и образования "провала" в направлении помехи
imagesc(angles_az*180/pi,angles_el*180/pi,pattern) % входной сигнал
title("ДНА")
ylabel("Логарифм амлитуды, норм.")
xlabel("Углы, град")
figure
imagesc(angles_az*180/pi,angles_el*180/pi,pattern_filtered) % отфильтрованный сигнал
title("ДНА после АПФ")
ylabel("Логарифм амлитуды, норм.")
xlabel("Углы, град")
%Вывод весовых коэффициентов
figure
imagesc(reshape(abs(weights),num_elem_az,num_elem_el))
%Вывод сечения в азимутальной плоскости
figure;
plot(angles_az*180/pi,pow2db(pattern(:,90)));hold on % Входной сигнал
plot(angles_az*180/pi,pow2db(pattern_filtered(:,90))); % Входной сигнал с ВК
ylim([-10,50])
grid on
title("ДНА в азимутальной плоскости")
ylabel("Логарифм амлитуды, норм.")
xlabel("Углы, град")
xline(az_signal/pi*180,'g') % задание на направление сканирования
legend ('ДНА без АПФ','ДНА после АПФ с функцией kaiser','Направление на сканирование','Направление на помеху','Направление на помеху')
for i=1:numel(JAM_DIR(:,1)) % задание направления на помеху
xline(JAM_DIR(i,1)/pi*180,'r')
end
% Вывод сечения в угломестной плоскости
figure;
plot(angles_el*180/pi,pow2db(pattern(90,:)));hold on % Входной сигнал
plot(angles_el*180/pi,pow2db(pattern_filtered(90,:))); % Входной сигнал с ВК
ylim([-10,50])
grid on
title("ДНА в угломестной плоскости")
ylabel("Логарифм амлитуды, норм.")
xlabel("Углы, град")
xline(el_signal/pi*180,'g') % задание на направление сканирования
legend ('ДНА без АПФ','ДНА после АПФ с функцией kaiser','Направление на сканирование')
for i=1:numel(JAM_DIR(:,2)) % задание направления на помеху
    xline(JAM_DIR(i,2)/pi*180,'r')
end
%% Вектор фазового набега
function[elem_phase_scan]=elem_phase_scan_func(angles,scan_direction,elem_pos)
elem_phase_scan = exp((sin(abs(angles-scan_direction))*elem_pos-...
fix(sin(abs(angles-scan_direction))*elem_pos))*2*pi*(-1i));
end
%% Расчет ДНА
% Расчёт разницы хода волн на элементах при электронном сканировании,
% разница хода волн - расстояние, нормированное к длине волны
function[pattern_calc]=pattern_calc_func(angles,weights,scan_direction,num_elem_az,num_elem_el,elem_pos,k_u)
num_elem = num_elem_az*num_elem_el;
elem_pattern = k_u*num_elem*prod(sinc(angles*pi/180));
elem_phase_scan = exp((sin(abs(angles-scan_direction))-...
fix(sin(abs(angles-scan_direction))))*elem_pos*2*pi*(1i));
pattern_calc = sum(1/num_elem.*...
elem_pattern*elem_phase_scan*weights);
end