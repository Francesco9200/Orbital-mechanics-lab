function [r_n, v_n] = propagator(t, dt, r_0, v_0, mu)
time = 0 : dt : t;
N = length(time);
r_n = [r_0];
v_n = [v_0];
for i = 1:N-1
   a1 = [v_n(:, i), -(mu/((norm(r_n(:, i)))^3))*r_n(:, i)];
   a2 = [v_n(:, i) + (dt/2)*a1(:, 2), -(mu/((norm((r_n(:, i) + a1(:, 1)*(dt/2))))^3))*(r_n(:, i) + a1(:, 1)*(dt/2))];
   a3 = [v_n(:, i) + (dt/2)*a2(:, 2), -(mu/((norm((r_n(:, i) + a2(:, 1)*(dt/2))))^3))*(r_n(:, i) + a2(:, 1)*(dt/2))];
   a4 = [v_n(:, i) + dt*a3(:, 2), -(mu/((norm((r_n(:, i) + a3(: ,1)*dt))))^3)*(r_n(:, i) + a3(:, 1)*dt)];
   r_n(:, i+1) = r_n(:, i) + dt/6*(a1(:, 1) + 2*a2(:, 1) + 2*a3(:, 1) + a4(:, 1));
   v_n(:, i+1) = v_n(:, i) + dt/6*(a1(:, 2) + 2*a2(:, 2) + 2*a3(:, 2) + a4(:, 2));
end
end
   