function tf = recSTFT(x, t, f, B)
  T = numel(t);
  F = numel(f);
  dt = t(2) - t(1);
  Q = B / dt;
  padx = [zeros(1, Q), x, zeros(1, Q)];
  padt = [zeros(1, Q), t, zeros(1, Q)];
  tf = zeros(T, F);
  
  % implement as cumulative sum can solve the accumulated error problem
  for m = 1: F
    x1 = padx .* exp(-2j*pi*padt*f(m)) * dt;
    cx1 = cumsum(x1);
    tf(:, m+1) = [cx1(2*Q+1), cx1(2*Q+2:2*Q+T) - cx1(1:T-1)];
  end
end