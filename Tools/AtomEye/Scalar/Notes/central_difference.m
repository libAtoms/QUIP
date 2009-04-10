%%%%%%%%%%%%%%%%%%%%%%%%
% central_difference.m %
%%%%%%%%%%%%%%%%%%%%%%%%

fprintf (1, 'central difference using f[-n],...f[n]:\n');
n = input('n = ');
fprintf (1, '%d-point central difference schemes:\n', 2*n+1);
                       
% odd derivative using combination (f[1]-f[-1]), .., (f[n]-f[-n]) %
% f^d (d even) derivatives are all canceled out %
% f^d (d odd) derivatives leaves f^d(2h^d/d!)n^d, .., f^d(2h^d/d!)1^d %
% we have n-1 degrees of freedom, therefore can choose to cancel the %
% the lowest n odd derivatives except the one we want to keep, dp %

% even derivative using combination %
% [(f[1]-f[0])+(f[-1]-f[0])], .., [(f[n]-f[0])+(f[-n]-f[0])] %
% f^d (d odd) derivatives are all canceled out %
% f^d (d even) derivatives leaves f^d(2h^d/d!)n^d, .., f^d(2h^d/d!)1^d %
% we have n-1 degrees of freedom, therefore can choose to cancel the %
% the lowest n odd derivatives except the one we want to keep, dp %

A = []; B = [];
for k = 1 : n,
  A(k,:) = (1:n) .^ (2*k-1);
  B(k,:) = (1:n) .^ (2*k);
end;

for k = 1 : n,
  b = zeros(n,1);
  b(k) = det(A);
  c = round(inv(A)*b);
  f = b(k);
  for l = 1 : n;
    f = gcd(f,c(l));
  end;
  c=c/f; b(k)=b(k)/f;
  d = 2*k - 1;
  fprintf (1, 'f^%d[0] = {', d);
  for l = 1 : n;
    fprintf (1, ' (%d)*(f[%d]-f[%d])', c(l), l, -l);
    if (l~=n) fprintf (1, ' +'); end;
  end;
  fprintf (1, ' } / %d / h^%d + O(h^%d)\n', ...
           b(k)/factorial(d)*2, d, 2*n-2*k+2);
  b(k) = det(B);
  c = round(inv(B)*b);
  f = b(k);
  for l = 1 : n;
    f = gcd(f,c(l));
  end;
  c=c/f; b(k)=b(k)/f;
  d = 2*k;
  fprintf (1, 'f^%d[0] = {', d);
  for l = 1 : n;
    fprintf (1, ' (%d)*((f[%d]-f[0])+(f[%d]-f[0]))', c(l), l, -l);
    if (l~=n) fprintf (1, ' +'); end;
  end;
  fprintf (1, ' } / %d / h^%d + O(h^%d)\n', ...
           b(k)/factorial(d)*2, d, 2*n-2*k+2);
end;
