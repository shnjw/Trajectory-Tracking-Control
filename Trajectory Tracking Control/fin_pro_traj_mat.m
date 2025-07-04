clc;
clear all;
close all;

%% 식별한 모델 매개변수
% >> 아래 매개변수 8개를 직접 채워 넣고 실험 진행할 것

hat_a11 = 2.104; 
hat_a12 = 1.992;
hat_a21 = 0.613;
hat_a22 = 1.302;

hat_b11 = 1.051;
hat_b22 = 1.521;

hat_f1 = 0.900;
hat_f2 = 2.159;

hat_A = [hat_a11 hat_a12; hat_a21 hat_a22];
hat_B = [hat_b11 0; 0 hat_b22];

%% 제어기 매개변수

Bi = inv(hat_B);
Kp = 50;
Kd = 50;
Ki = 10;
lamda = 100;
k1 = 100;
k2 = 100;


%% 
out = sim("fin_pro_traj_sim.slx");

% 샘플링 시간과 유효 데이터 수
dt = 0.005;
t = out.e1.Time;  % [4001×1]
n = round(20 / dt);  % 20초 이내 길이 선택

% 데이터 1차원 벡터로 변환 (squeeze)
e1 = squeeze(out.e1.Data);  % [4001×1]
e2 = squeeze(out.e2.Data);  % 동일하게 적용

% MAE 계산
mae1 = mean(abs(e1));
mae2 = mean(abs(e2));

fprintf('MAE1 = %.6f\n', mae1);
fprintf('MAE2 = %.6f\n', mae2);








