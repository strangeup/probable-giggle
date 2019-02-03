% LOAD the trace
el_area=load("areas.dat")
l2_norm=load("errors.dat")
% Typical 'size' of elements is sqrt(element_area)
% Get the gradient and intercept
m=polyfit(0.5*log10(el_area),log10(l2_norm),1)
% Brewer blue
linestyle = linspecer(1,"qualitative")
% Begin a figure
figure;
hold on
xlabel 'log_{10}h_n'
ylabel 'log_{10} l^2'
% Plot the line and the fit
plot(0.5*log10(el_area),log10(l2_norm),'.','markerfacecolor',linestyle)
plot(0.5*log10(el_area),m(1)*0.5*log10(el_area)+m(2),'-','color',linestyle)
