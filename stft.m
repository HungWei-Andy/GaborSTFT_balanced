function [tf, elapsedTime] = stft(x, w, deltat, deltaf, basef, F, type)
% This function calls "naive_stft", "fft_stft", "rec_stft"
% The default frequence interval is [-F_half*deltaf : (F-F_half-1)*deltaf]
%
% x: sampled time signal
% w: sampled window 
% deltat: sampling interval on time, it's assumed that deltat = deltatau
% deltaf: sampling interval on frequency
% F: number of sampling points on frequency. That is, the sampling  points
%    are 0, deltaf, deltaf * 2, ... , deltaf * (F-1).
% basef: base frequency. the real sampling points are basef, basef +
%        deltaf, basef + deltaf*2, ..., basef + deltaf * (F-1).

  % check
  if( ~isvector(x) || ~isvector(w) ) 
    throw MException('MyExcpetion: type error', 'x and w should be vectors');
  end
  if( ~isscalar(deltat) || ~isscalar(deltaf) || ~isscalar(F) )
    throw MException('MyExcpetion: type error', 'deltat, deltaf, F should be scalars');
  end
  if( ~isnumeric(x) || ~isnumeric(w) || ~isnumeric(deltat) || ~isnumeric(deltaf) || ~isnumeric(F) )
    throw MException('MyException: numeric error', 'x, w, deltat, deltaf, F should be numeric';
  end
  
  % w should contain odd number of points; otherwise, pad a zero
  if( mod(numel(w), 2) == 0 )
    w = [0 w];
  end
  
  % compute
  if( strcmp(type, 'naive') )
    if(nargout == 2)
      [tf, elapsedTime] = naive_stft(x, w, deltat, deltaf, basef, F);
    elseif(nargout == 1)
      tf = naive_stft(x, w, deltat, deltaf, basef, F);
    end
  elseif( strcmp(type, 'fft') )
    if(nargout == 2)
      [tf, elapsedTime] = fft_stft(x, w, deltat, deltaf, basef, F);
    elseif(nargout == 1)
      tf = fft_stft(x, w, deltat, deltaf, basef, F);
    end
  elseif( strcmp(type, 'recursive') )
    if(nargout == 2)
      [tf, elapsedTime] = rec_stft(x, numel(w), deltat, deltaf, basef, F);
    elseif(nargout == 1)
      tf = rec_stft(x, numel(w), deltat, deltaf, basef, F);
    end
  else
    throw MException('MyException: method not exist', ['no such method type: ', type]);
  end
end

% -------------------- implementation functions ---------------------------

function [tf, elapsedTime] = naive_stft(x, w, deltat, deltaf, basef, F)
  % pad x with zeros
  T = numel(x);
  Q = (numel(w) - 1) / 2;
  padx = [zeros(1, Q), x, zeros(1, Q)];
  tf = zeros(T, F);
  
  % start timer
  if(nargout == 2)
    tic();
  end
  
  % compute
  for n = 0 : T-1
    for m = 0 : F-1
      tf(n+1, m+1) = sum(fliplr(w) .* padx(n+1: n+2*Q+1) .* exp(-1i*2*pi*(n-Q:n+Q)*deltat*(m*deltaf + basef)) * deltat);
    end
  end
  
  % end timer
  if(nargout == 2)
    elapsedTime = toc();
  end
end

function [tf, elapsedTime] = fft_stft(x, w, deltat, deltaf, basef, F)
  % check constraints
  % 1. deltat * deltaf ~= integer
  % 2. N >= 2*Q+1
  % 3. N >= F
  T = numel(x);
  Q = (numel(w) - 1) / 2;
  N = round(1/(deltat*deltaf));
  if( abs(1/(deltat*deltaf) - N) > 1e-3 )
    throw MException('MyException: fft constraints violated', 'deltat * deltaf ~= 1/N');
  end
  if( N < 2*Q+1 )
    throw MException('MyException: fft constraints violated', 'N >= 2Q+1');
  end
  if( N < F)
    throw MException('MyException: fft constraintes violated', 'N >= F');
  end
  
  %arguments
  padx = [zeros(1, Q), x, zeros(1, Q)];
  tf = zeros(T, F);
  
  % start timer
  if(nargout == 2)
    tic();
  end
  
  % compute
  for n = 0 : T-1
    x1 = [fliplr(w) .* padx(n+1: n+2*Q+1) .* exp(-1i*2*pi*(n-Q:n+Q)*deltat*basef), zeros(1, N-2*Q-1)];
    fftx1 = fft(x1);
    tf(n+1, :) = fftx1(1: F) .* exp(1i*2*pi*(Q-n)*(0:F-1)) * deltat;
  end
  
  % end timer
  if(nargout == 2)
    elapsedTime = toc();
  end
end

function [tf, elapsedTime] = rec_stft(x, wlen, deltat, deltaf, basef, F)
  T = numel(x);
  Q = (wlen - 1) / 2;
  padx = [zeros(1, Q), x, zeros(1, Q)];
  tf = zeros(T, F);
  
  % start timer
  if(nargout == 2)
    tic();
  end
  
  for m = 0: F-1
    % implemented as difference of cumlative sum can solve the accumulated
    % error problem
    x1 = padx .* exp(-1i*2*pi*(-Q:T-1+Q)*deltat*(m*deltaf+basef)) * deltat;
    cx1 = cumsum(x1);
    tf(:, m+1) = [cx1(2*Q+1), cx1(2*Q+2:2*Q+T) - cx1(1:T-1)];
  end
  
  % end timer
  if(nargout == 2)
    elapsedTime = toc();
  end
end