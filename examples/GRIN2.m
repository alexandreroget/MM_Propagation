filename = 'input_GRIN2';

%% Fiber parameters

n_modes = 2;

fiber_length = 16; % in m

max_dispersion_order = 2;

dispersion_coefficients = zeros(n_modes,max_dispersion_order+1);
dispersion_coefficients(1,3) = 18.938e-3;
dispersion_coefficients(2,3) = 18.866e-3;

alpha = 6.4341e12;
coupling_coefficients = zeros(2,2,2,2);

coupling_coefficients(1,1,1,1) = alpha;
coupling_coefficients(1,1,1,2) = alpha/3;
coupling_coefficients(1,1,2,2) = alpha/2;
coupling_coefficients(1,2,2,2) = -alpha/3;

coupling_coefficients(2,1,1,1) = -alpha/3;
coupling_coefficients(2,2,1,1) = alpha/5;
coupling_coefficients(2,2,2,1) = alpha/3;
coupling_coefficients(2,2,2,2) = 0.8*alpha;

%% Simulation parameters

n_steps = 131072;

order_RK = 4;

speed_of_light = 2.99792458e-4; % m/ps
lambda = 1030e-9; % m, the center wavelength
f0 = speed_of_light/lambda;
w0 = 2*pi*f0; % angular frequency (THz)
n2 = 2.3*10^-20; % m^2 W^-1
nonlinearity_const = n2*w0/speed_of_light; % W^-1 m

raman_proportion = 0;
raman_parameters = [12e-3 , 32e-3];

%% Initial conditions

N = 2^18; % Number of time/frequency points
time_window = 160; % ps
pulse_width = sqrt(dispersion_coefficients(1,3))/(alpha*nonlinearity_const);

% t0 = -time_window/10;
t0 = 0;

q = alpha*nonlinearity_const;
c = 1 * sqrt(2/dispersion_coefficients(1,3));

dt = time_window/N;
t = linspace(-N/2,N/2-1,N).*dt;

initial_fields = zeros(N,2);
initial_fields(:,1) = (0.5 * sqrt(q/2)) .* sech(q/4 * (t-t0)/sqrt(dispersion_coefficients(1,3)/2)) .* exp((1i * c) .* (((t-t0)/sqrt(dispersion_coefficients(1,3)/2))./2));
initial_fields(:,2) = (0.5 * sqrt(0.8*q/2)) .* sech(0.8*q/4 * (t+t0)/sqrt(dispersion_coefficients(1,3)/2)) .* exp((-1i * c) .* (((t+t0)/sqrt(dispersion_coefficients(1,3)/2))./2));

%% Save parameters

save(filename,'n_modes','fiber_length','max_dispersion_order','dispersion_coefficients','coupling_coefficients','N','time_window','pulse_width','initial_fields','n_steps','order_RK','nonlinearity_const','raman_proportion','raman_parameters');