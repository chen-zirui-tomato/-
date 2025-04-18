% 定义 Adams-Bashforth 的 ρ 和 σ 多项式（p=1 到 5）
orders = 1:5;
theta = linspace(0, 2*pi, 1000);
z = exp(1i*theta);

figure;
hold on;
for p = orders
    % 根据阶数 p 定义 ρ(z) 和 σ(z) 的系数
    switch p
        case 1
            rho = [1, -1];       % ρ(z) = z - 1
            sigma = [0, 1];      % σ(z) = 1
        case 2
            rho = [1, -1, 0];    % ρ(z) = z^2 - z
            sigma = [0, 3/2, -1/2];
        case 3
            rho = [1, -1, 0, 0];
            sigma = [0, 23/12, -16/12, 5/12];
        case 4
            rho = [1, -1, 0, 0, 0];
            sigma = [0, 55/24, -59/24, 37/24, -9/24];
        case 5
            rho = [1, -1, 0, 0, 0, 0];
            sigma = [0, 1901/720, -2774/720, 2616/720, -1274/720, 251/720];
    end
    
    % 计算 γ(θ) = ρ(z)/σ(z)
    gamma = polyval(rho, z) ./ polyval(sigma, z);
    
    % 绘制根轨迹曲线
    plot(real(gamma), imag(gamma), 'LineWidth', 1.5, 'DisplayName', ['p=', num2str(p)]);
end

xlabel('Re(κ)'); ylabel('Im(κ)');
title('Adams-Bashforth 绝对稳定域 (p=1 到 5)');
legend('Location', 'best');
grid on;