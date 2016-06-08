figure(1);
plot(segment,freq,'b',segment,-1*freq,'b');
xlabel('Index of Forced Oscillator');
ylabel('Frequency Difference \omega_f-\omega_0 (rad/s)');
title('Entrainment Ranges for Non-Uniform Coupling Asymmetry');
screen_size = get(0, 'ScreenSize');
set(1, 'Position', [0 0 0.75*screen_size(3) 0.5*screen_size(4) ] );
export_fig('EntrainmentRangeNonUniform','-pdf','-nocrop');