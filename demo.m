deltat = 1e-2;
B = 0.6;
f = -5 : 0.05 : 5;
t = [0: deltat: 15-deltat];

x1 = cos(2*pi*(0: deltat: 5-deltat));
x2 = cos(7*pi*(5: deltat: 10-deltat));
x3 = cos(4*pi*(10: deltat: 15-deltat)); 
x = [x1, x2, x3];

T = numel(t);

% STFT
tic();
tf = recSTFT(x, t, f, B);
elapsedTime = toc();
disp(['elapsed Time: ' elapsedTime]);
  
% find norm and convert into 500x500 image 
tf_norm = abs(tf);
img = imresize(tf_norm', [500, 500]);

% save into images
h = figure;
imshow(img, ...
      'XData', t, ...
      'YData', f);
set(gca,'YDir','normal');
title(['Chirp Function' ', elapsed time = ' num2str(elapsedTime) 'sec, matrix size (time freq)= ' num2str(size(tf_norm))]);
xlabel('time (sec)');
ylabel('frequency (Hz)');
axis on;
  
% save the figure as ".fig" and ".jpg"
savefig(h, 'output.fig');
saveas(h, 'output.jpg');
