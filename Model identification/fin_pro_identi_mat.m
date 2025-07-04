clc;
clear all;
close all;




%% 설정
dim = 8;  % a1, a2, a3, b1, b2, b3
iter_max = 80;
lambda = 40; % 자식 개체 수
mu = floor(lambda / 2); % 부모 개체 수
lu = [0; 10]; % lower & upper bound

% 초기화
xmean = ones(1, dim);  % 초기 평균
sigma = 0.8;               % 초기 표준편차
weights = log(mu + 0.5) - log(1:mu); % weight: 우수 해 가중치 
weights = weights / sum(weights); 
mueff = sum(weights)^2 / sum(weights.^2); % mueff: 실질적 부모 개체 수

% CMA 파라미터
cc = (4 + mueff/dim) / (dim + 4 + 2*mueff/dim);% 공분산 업데이트 학습률 p_c
cs = (mueff + 2) / (dim + mueff + 5); % 스텝사이즈 학습률 p_sigma
c1 = 2 / ((dim + 1.3)^2 + mueff); % Rank-1 업데이트 비율 
cmu = min(1 - c1, 2 * (mueff - 2 + 1/mueff) / ((dim + 2)^2 + mueff));% Rank-u 업데이트 비율
damps = 1 + 2 * max(0, sqrt((mueff - 1)/(dim + 1)) - 1) + cs;

pc = zeros(1, dim);
ps = zeros(1, dim);
B = eye(dim);
D = ones(dim, 1);
C = B * diag(D.^2) * B';
eigeneval = 0;
chiN = dim^0.5 * (1 - 1/(4*dim) + 1/(21*dim^2));

cost_REC = zeros(1, iter_max);

% 데이터 로딩
%% data

load data1.mat % data 1~10까지 있음

x1 = out.x1;
x2 = out.x2;

u1.time = out.tout;
u1.signals.dimensions = 1;
u1.signals.values = out.u1;

u2.time = out.tout;
u2.signals.dimensions = 1;
u2.signals.values = out.u2;


for gen = 1:iter_max
    % 샘플링
    arz = randn(lambda, dim);             % 정규분포
    ary = arz * (B * diag(D))';           % C = B D^2 B'
    arx = xmean + sigma * ary;            % 변형된 샘플
    arx = min(max(arx, lu(1,:)), lu(2,:)); % 경계 제약

    % 평가
    fitness = zeros(lambda, 1);
    for k = 1:lambda
        hat_A = ([arx(k,1),arx(k,2); arx(k,3), arx(k,4)]);
        hat_B = ([arx(k,5), 0; 0, arx(k,6)]);
        hat_f1 = arx(k,7);
        hat_f2 = arx(k,8);
        from_sim = sim("fin_pro_identi_sim.slx");
        hat_x1 = from_sim.hat_x1;
        hat_x2 = from_sim.hat_x2;
        fitness(k) = sum(sum(abs(x1 - hat_x1)) + sum(abs(x2 - hat_x2)));
    end

    % 정렬 및 선택
    [~, idx] = sort(fitness);
    xold = xmean;
    xmean = weights * arx(idx(1:mu), :);    % 새로운 평균

    % 업데이트
    y = (xmean - xold) / sigma;
    ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * (y / (B * diag(D))');
    hsig = norm(ps) / sqrt(1 - (1 - cs)^(2 * gen)) / chiN < 1.4 + 2 / (dim + 1);
    pc = (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * y;

    % 공분산 행렬 업데이트
    artmp = (arx(idx(1:mu), :) - xold) / sigma;
    C = (1 - c1 - cmu) * C + ...
        c1 * (pc' * pc + (1 - hsig) * cc * (2 - cc) * C) + ...
        cmu * artmp' * diag(weights) * artmp;

    sigma = sigma * exp((cs / damps) * (norm(ps)/chiN - 1));

    % 고유값 분해
    if mod(gen, 1) == 0
        C = triu(C) + triu(C, 1)'; 
        [B, D_diag] = eig(C); 
        D = sqrt(diag(D_diag)); 
    end

    cost_REC(gen) = fitness(idx(1));
    disp(['Iter ' num2str(gen) ', Best Cost = ' num2str(cost_REC(gen))])
end

%% 결과 시각화

figure();
plot(cost_REC, 'LineWidth', 2);
xlabel("Iteration"); ylabel("Cost");
title("CMA-ES 식별 Cost 수렴");
grid on;

%% 결과 확인
params = xmean;
hat_A = ([params(1), params(2); params(3), params(4)]);
hat_B = ([params(5), 0; 0, params(6)]);
hat_f1 = params(7);
hat_f2 = params(8);
from_sim = sim("fin_pro_identi_sim.slx");
hat_x1 = from_sim.hat_x1;
hat_x2 = from_sim.hat_x2;


figure()
hold on
title('$$x1$$', 'Interpreter','latex')
plot(out.tout, x1)
plot(from_sim.tout, from_sim.hat_x1)
legend('$$x1$$','$$hatx1$$', 'Interpreter','latex', 'Box','off')
xlabel("time [s]")
ylabel("angle [rad]")
ax = gca; ax.FontSize = 15;

figure()
hold on
title('$$x2$$', 'Interpreter','latex')
plot(out.tout, x2)
plot(from_sim.tout, from_sim.hat_x2)
legend('$$x2$$','$$hatx2$$', 'Interpreter','latex', 'Box','off')
xlabel("time [s]")
ylabel("angle [rad]")
ax = gca; ax.FontSize = 15;