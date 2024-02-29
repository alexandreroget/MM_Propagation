%%%%%%%%%%%%%%%%%%%%%%%% Fiber parameters %%%%%%%%%%%%%%%%%%%%%%%%%

n_modes = 2;

fiber_length = 16; % in m

order_RK = 4;

dispersion_coefficients = zeros(n_modes,3);
dispersion_coefficients(1,3) = 18.938e-3;
dispersion_coefficients(2,3) = 18.866e-3;

alpha = 6.4341e12;
coupling_coefficients = zeros(2,2,2,2);

coupling_coefficients(1,1,1,1) = alpha;
coupling_coefficients(1,1,1,2) = alpha/3;
coupling_coefficients(1,1,2,2) = alpha/2;
% coupling_coefficients(1,2,2,2) = -alpha/3;

% coupling_coefficients(2,1,1,1) = -alpha/3;
coupling_coefficients(2,2,1,1) = alpha/5;
coupling_coefficients(2,2,2,1) = alpha/3;
coupling_coefficients(2,2,2,2) = 0.8*alpha;

%%%%%%%%%%%%%%%%%%%%%%% Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%

speed_of_light = 2.99792458e-4; % m/ps
lambda = 1030e-9; % m, the center wavelength
f0 = speed_of_light/lambda;
w0 = 2*pi*f0; % angular frequency (THz)
n2 = 2.3*10^-20; % m^2 W^-1
gamma = n2*w0/speed_of_light; % W^-1 m

N = 2^15; % Number of time/frequency points
time_window = 160; % ps
pulse_width = (sqrt(dispersion_coefficients(1,3))/(alpha*gamma));

t0 = -time_window/10;

q = alpha*gamma;
c = 1 * sqrt(2/dispersion_coefficients(1,3));

dt = time_window/N;
t = linspace(-N/2,N/2-1,N).*dt;

initial_fields = zeros(2,N);
initial_fields(1,:) = (0.5 * sqrt(q/2)) .* sech(q/4 * (t-t0)/sqrt(dispersion_coefficients(1,3)/2)) .* exp((1i * c) .* (((t-t0)/sqrt(dispersion_coefficients(1,3)/2))./2));
initial_fields(2,:) = (0.5 * sqrt(0.8*q/2)) .* sech(0.8*q/4 * (t+t0)/sqrt(dispersion_coefficients(1,3)/2)) .* exp((-1i * c) .* (((t+t0)/sqrt(dispersion_coefficients(1,3)/2))./2));

%%%%%%%%%%%%%%%%%%%%%% Simulation parameters %%%%%%%%%%%%%%%%%%%%%%

n_steps = 1024;

nonlinearity_const = gamma;

raman_proportion = 0.18;

t1 = 12.2e-3; % raman parameter t1 [ps]
t2 = 32e-3; % raman parameter t2 [ps]
t_shifted = linspace(0,N-1,N).*dt;
raman_response = ((t1^2+t2^2)/(t1*t2^2)).*exp(-t_shifted/t2).*sin(t_shifted/t1);
raman_response = transpose(raman_response);

%%%%%%%%%%%%%%%%%%%%%%%%% Save parameters %%%%%%%%%%%%%%%%%%%%%%%%%

input.n_modes = n_modes;
input.fiber_length = fiber_length;
input.dispersion_coefficients = dispersion_coefficients;
input.coupling_coefficients = coupling_coefficients;
input.initial_fields = initial_fields;
input.time_window = time_window;
input.pulse_width = pulse_width;
input.n_steps = n_steps;
input.method_order = 4;
input.nonlinearity_const = nonlinearity_const;
input.savename = 'GRIN2';
input.n_save = 4;
input.raman_proportion = raman_proportion;
input.raman_response = raman_response;

output = mm_propagation(input);