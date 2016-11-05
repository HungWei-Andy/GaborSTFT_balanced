function show_stft(x, w, deltat, deltaf, basef, F, type, msg)
% This function calls stft to generate the time-frequency spectrum and 
% find its norm by calling abs() and display the norm as an  image.
  T = numel(x);
  
  % STFT
  [tf, elapsedTime] = stft(x, w, deltat, deltaf, basef, F, type);
  
  % find norm and convert into 500x500 image 
  tf_norm = abs(tf);
  img = imresize(tf_norm', [500, 500]);
  
  % save into images
  h = figure;
  set(h, 'Visible', 'off'); % don't display the figure
  imshow(img, ...
        'XData', 0: T/10*deltat: deltat * T, ...
        'YData', basef: F/10*deltaf: basef + F*deltaf);
  set(gca,'YDir','normal');
  title([msg ', elapsed time = ' num2str(elapsedTime) 'sec, matrix size (time freq)= ' num2str(size(tf_norm))]);
  xlabel('time (sec)');
  ylabel('frequency (Hz)');
  axis on;
  colorbar;
  
  % save the figure as ".fig" and ".jpg"
  savefig(h, ['res/' msg]);
  saveas(h, ['res/' msg '.jpg']);
end