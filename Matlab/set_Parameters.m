filename = 'input_MMPropagation';

%% Fiber parameters

n_modes = 4;

fiber_length = 100; % in m

max_dispersion_order = 2;

dispersion_coefficients = zeros(n_modes,max_dispersion_order+1);

coupling_coefficients = zeros(n_modes,n_modes,n_modes,n_modes);

%% Simulation parameters

n_steps = 1024;

order_RK = 4;

speed_of_light = 2.99792458e-4; % m/ps
lambda = 1030e-9; % m, the center wavelength
f0 = speed_of_light/lambda;
w0 = 2*pi*f0; % angular frequency (THz)
n2 = 2.3*10^-20; % m^2 W^-1
nonlinearity_const = n2*w0/speed_of_light; % W^-1 m

raman_proportion = 0.18;
raman_parameters = [12e-3 , 32e-3];

%% Initial conditions

N = 2^20; % Number of time/frequency points
time_window = 160; % ps
pulse_width = 2;

initial_fields = zeros(N,n_modes);

%% Save parameters

save(filename,'n_modes','fiber_length','max_dispersion_order','dispersion_coefficients','coupling_coefficients','N','time_window','pulse_width','initial_fields','n_steps','order_RK','nonlinearity_const','raman_proportion','raman_parameters');