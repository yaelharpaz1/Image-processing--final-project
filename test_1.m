function [t,u]=test_1()

%% Synthetic test image:
C = 0.5.*ones(256,256);
C(101:156,21:235) = 0.3;
C(21:235,101:156) = 0.3;
I = mat2gray(C,[1 0]);
figure
imshow(I)
title ('Test image')

%I=imread('lena_color.tiff');
%I=rgb2gray(I);
%figure
%imshow(I)
%title ('Test image')

ref = I;

% Adding noise:
u = imnoise(I,'gaussian');
figure
imshow(u)                % Adding gaussian white noise with var of 0.01
av_GL = mean(u(:));
title (sprintf('t = 0  Average Gray level = %.3g',av_GL)) 
[peaksnr, snr] = psnr(u, ref);
xlabel (sprintf('PSNR = %.3g  SNR = %.3g',peaksnr,snr));          

%% Pixel size:
h_1 = 1;
h_2 = 1;

%% The diffusivity function- Type 2:
k = 2.5;                % k>1
lambda = 2;             % lambda>0
g = @(s) 2*exp(-((k^2*log(2))/(k^2-1))*(s/lambda^2)) - ...
        exp(-(log(2)/(k^2-1))*(s/lambda^2));

R = 1;                  % grey values interval length.
c_1 = g(0);             % the diffusivity extremum.
[~,c_2] = fminbnd(g,0,10^4);

syms wR
eqn = g(wR) == -c_2;
w = solve(eqn,wR)/R;    % the stabilisation range constant.

%% Max-Min Principle- step size:
t_min = (w^2*h_1^4*h_2^4)/(2*c_1*(h_1^2 + h_2^2)*(w^2*h_1^2*h_2^2 + h_1^2 + h_2^2));
t_max = 1/(2*c_1*(1/h_1^2 + 1/h_2^2));

%%
t = 0;                  % starting time.
t_final = 100;          % stopping time.

while t<10
    [u,tau] = FAB_GlobalStepSize(u,t_max,t_min,h_1,h_2,g);
    t = t + tau;
    disp(t)
end 

figure
imshow(u)
av_GL = mean(u(:));
title (sprintf('t = %.3g  Average Gray level = %.3g',t,av_GL))
[peaksnr, snr] = psnr(u, ref);
xlabel (sprintf('PSNR = %.3g  SNR = %.3g',peaksnr,snr));

while t<20
    [u,tau] = FAB_GlobalStepSize(u,t_max,t_min,h_1,h_2,g);
    t = t + tau;
    disp(t)
end 

figure
imshow(u)
av_GL = mean(u(:));
title (sprintf('t = %.3g  Average Gray level = %.3g',t,av_GL))
[peaksnr, snr] = psnr(u, ref);
xlabel (sprintf('PSNR = %.3g  SNR = %.3g',peaksnr,snr));

while t<40
    [u,tau] = FAB_GlobalStepSize(u,t_max,t_min,h_1,h_2,g);
    t = t + tau;
    disp(t)
end 

figure
imshow(u)
av_GL = mean(u(:));
title (sprintf('t = %.3g  Average Gray level = %.3g',t,av_GL))
[peaksnr, snr] = psnr(u, ref);
xlabel (sprintf('PSNR = %.3g  SNR = %.3g',peaksnr,snr));

while t<60
    [u,tau] = FAB_GlobalStepSize(u,t_max,t_min,h_1,h_2,g);
    t = t + tau;
    disp(t)
end 

figure
imshow(u)
av_GL = mean(u(:));
title (sprintf('t = %.3g  Average Gray level = %.3g',t,av_GL))
[peaksnr, snr] = psnr(u, ref);
xlabel (sprintf('PSNR = %.3g  SNR = %.3g',peaksnr,snr));

while t<80
    [u,tau] = FAB_GlobalStepSize(u,t_max,t_min,h_1,h_2,g);
    t = t + tau;
    disp(t)
end 

figure
imshow(u)
av_GL = mean(u(:));
title (sprintf('t = %.3g  Average Gray level = %.3g',t,av_GL))
[peaksnr, snr] = psnr(u, ref);
xlabel (sprintf('PSNR = %.3g  SNR = %.3g',peaksnr,snr));

while t<t_final
    [u,tau] = FAB_GlobalStepSize(u,t_max,t_min,h_1,h_2,g);
    t = t + tau;
    disp(t)
end 

figure
imshow(u)
av_GL = mean(u(:));
title (sprintf('t = %.3g  Average Gray level = %.3g',t,av_GL))
[peaksnr, snr] = psnr(u, ref);
xlabel (sprintf('PSNR = %.3g  SNR = %.3g',peaksnr,snr));
