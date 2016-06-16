% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20)

figure(1)
load Nonuniform.mat
plot(segment,freq,'b',segment,-1*freq,'b', 'LineWidth', 3)
title('Nonuniform Coupling Asymmetry','FontSize', 20)
A = legend('Entrainment Threshold')
set(A,'Interpreter','latex')
set(gca,'fontsize',20)
ylabel('Frequency Difference wf-w0 (rad/s)','FontSize', 20)
xlabel('Index of Forced Oscillator m','FontSize', 20)

figure(2)
load Uniform.mat
plot(segment,freq,'b',segment,-1*freq,'b', 'LineWidth', 3)
title('Uniform Coupling Asymmetry','FontSize', 20)
B = legend('Entrainment Threshold')
set(B,'Interpreter','latex')
set(gca,'fontsize',20)
ylabel('Frequency Difference wf-w0 (rad/s)','FontSize', 20)
xlabel('Index of Forced Oscillator m','FontSize', 20)