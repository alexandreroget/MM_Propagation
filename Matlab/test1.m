filename = 'test_RK4_grand';

%% Fiber parameters

n_modes = 1;

fiber_length = 100; % in m

max_dispersion_order = 2;

dispersion_coefficients = ones(n_modes,max_dispersion_order+1);

coupling_coefficients = 1;

%% Simulation parameters

n_steps = 2048;

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

N = 2^15; % Number of time/frequency points
time_window = 160; % ps
pulse_width = 2;

t_final = time_window;
t = linspace(0,t_final,N);
initial_fields = transpose(exp(1i.*(t./t_final)));

%% Save parameters

save(filename,'n_modes','fiber_length','max_dispersion_order','dispersion_coefficients','coupling_coefficients','N','time_window','pulse_width','initial_fields','n_steps','order_RK','nonlinearity_const','raman_proportion','raman_parameters');