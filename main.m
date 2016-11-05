deltat = 1e-2;
deltaf = 0.05;
F = 200;

% examples in ppt
% window
gaborQ = 2;
gabor_w = exp(-pi*(-gaborQ: deltat: gaborQ).^2); % Gabor Transform
recQ = 0.5;
rec_w = ones(1, 2*recQ/deltat + 1); % Rectangular STFT

% signal
x1 = exp(-pi*(-10: deltat: 10).^2); % Gaussian Function
basef1 = -5;

x2_1 = cos(2*pi*(0: deltat: 4-deltat));
x2_2 = cos(6*pi*(4: deltat: 8-deltat));
x2_3 = cos(4*pi*(8: deltat: 12-deltat));  % Cosine Function
x2 = [x2_1, x2_2, x2_3];
basef2 = 0;

t = -9: deltat: 1;
x3 = [zeros(1, 10), exp(1i*t.^2/10 - 1i*3*t), zeros(1, 10)]; % Chirp Function
basef3 = -5;

t = -9: deltat: 9;
x4 = exp(1i*t.^2/2 + 1i*6*t) .* exp(-(t-4).^2/10); % Exponential Chirp Function
basef4 = -2;

% compute stft
ws = {gabor_w, rec_w};
xs = {x1, x2, x3, x4};
basefs = {basef1, basef2, basef3, basef4};
ws_name = {'Gabor window', 'Rect window'};
xs_name = {'Exp', 'Cosine', 'Chirp', 'Exp Chirp'};
types = {'naive', 'fft', 'recursive'};

for wi = 1: numel(ws)
  for xi = 1: numel(xs)
    for ti = 1: numel(types)
      % Gabor window, don't try recursive implementation
      if( strcmp(ws_name{wi}, 'Gabor window') && strcmp(types{ti}, 'recursive') )
        continue;
      end
      
      show_stft(xs{xi}, ws{wi}, deltat, deltaf, basefs{xi}, F, types{ti}, ...
          [xs_name{xi} ' (' ws_name{wi} ', ' types{ti} ')']);
    end
  end
end

