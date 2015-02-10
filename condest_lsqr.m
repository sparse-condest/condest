function [sigma_max,sigma_min,v_max,v_min,stats] = condest_lsqr(A, opts) % upper_bound, params.verify,trace_interval,x_orig)

fontsize = 14;
legend_fontsize = 12;
linewidth = 1.5;
markersize = 9;

if nargin < 1
    error 'No input matrix given'
end

[m, n] = size(A);
if (m < n)
    A = A';
    [m, n] = size(A);
end

params.type           = 'double';
params.shift          = 0;
params.trace_with_svd = false;
params.verify         = false;
params.trace          = false;
params.maxits         = 100000;
params.fixedits       = false;
params.extraits       = true;
params.restol         = 8*eps;
params.highcond_restol = 4*eps;
params.highcond_thresh = sqrt(eps);
% params.normdifftol    = sqrt(0.01) / sqrt(2*n);
params.c2 = 1e-3;

params.cond_thresh    = (1/eps)/64;
params.projtrace      = false;
OP=A;

if nargin >= 2
  if isfield(opts,'type');           params.type           = opts.type;           end
  if isfield(opts,'shift');          params.shift          = opts.shift;          end
  if isfield(opts,'trace_with_svd'); params.trace_with_svd = opts.trace_with_svd; end
  if isfield(opts,'verify');         params.verify         = opts.verify;         end
  if isfield(opts,'trace');          params.trace          = opts.trace;          end
  %if isfield(opts,'graph');          params.graph          = opts.graph;          end
    
  if isfield(opts,'maxits');         params.maxits         = opts.maxits;         end
  if isfield(opts,'fixedits');       params.fixedits       = opts.fixedits;       end
  if isfield(opts,'extraits');       params.extraits       = opts.extraits;       end
  if isfield(opts,'restol');         params.restol         = opts.restol;         end
  if isfield(opts,'c2');    params.c2    = opts.c2;       end
  if isfield(opts,'cond_thresh');        params.cond_thresh        = opts.cond_thresh;        end
  if isfield(opts,'highcond_thresh');        params.highcond_thresh        = opts.highcond_thresh;        end
  if isfield(opts,'highcond_restol');        params.highcond_restol        = opts.highcond_restol;        end
  if isfield(opts,'projtrace');      params.projtrace      = opts.projtrace;      end
  if isfield(opts,'precond');        OP=A*opts.precond;                           end
end

shift = numtype(params.shift);

A=numtype(A);
class(A)
if strcmp(class(A),'sym')
    disp('high-precision arithmetic');
end

if strcmp(params.type,'sym')
    params.verify = false;
    params.trace  = false;
end

vsteps = 5000; % maximum number of params.verify/trace steps
if params.trace
  %params.verify = true;
  trace_interval=1; 
end

% LSQR or RRLSQR?

explicit_r = false;
%explicit_r = true;

% estimate the 2-norm of A
%v = randn(n,1);
%v = v / norm2(v);
%for k=1:33 % based on the Klein-Lu STOC 1996 paper with epsilon=1e-1, delta=1e-12
%    u = A*v;
%    v = A'*u;
%    v = v / norm2(v);
%end
%v_max = v;
%sigma_max = norm2(A*v) / norm2(v);
[sigma_max,v_max] = power_rect_explicit(A);
sigma_max

A_orig = A;
if shift~=0
  A = [ shift*eye(n) ; A ];
  m = n+m;
end

% params.verify steps
if params.verify
  sA=svd(full(A)); 
  s_min = sA(n);
  normA = sA(1);
  fprintf(1,'norm2(A)=%.2e, power-method relative error=%.0e\n',double(sigma_max),double(abs(normA-sigma_max)/normA));
  
  U = zeros(m,vsteps);
  V = zeros(n,vsteps);
  %B = zeros(vsteps+1,vsteps);
  %R = zeros(vsteps,vsteps);
  P = zeros(m,vsteps);
end

% xxx for testing
if params.projtrace
[UA,SA,VA]=svd(A,0);
v_min_exact = VA(:,n);
if params.projtrace
    v_pt_min = VA(:,n);
    v_pt_mid = VA(:,round(n/1));
    v_pt_max = VA(:,1);
end
end
% xxx end

% Construct the shifted problem

sigma_min = sigma_max; % we set it to sigma_max to get a correct result on some matrices with cond=1
v_min = NaN*zeros(n,1);
normdiff_at_min = NaN;
iter_at_min = -1;
estimate_rising = false;

% Now try to solve min(A*xhat-b) and use xhat-x as an approximate null
% vector; we use LSQR

% LSQR Initialization
i = 1;
stop = false;
stop_res = NaN;
stop_normdiff = NaN;
stop_cond = NaN

% set up the consistent LSQR solve
% A representation for R
if nargin<5
  x_orig = numtype(randn(n,1));
end
normdifftol    = sqrt(2) * erfinv(params.c2) / norm2(x_orig);
x_orig = x_orig / norm2(x_orig);
b_c = A*x_orig;   % consistent
norm_b_c = norm2(b_c);
%maxabs_b = max(abs(b_c));
%pert = 8*eps*norm_b_c*(randn(m,1) > 0);
%size(pert)
%size(b_c)
%b_c = b_c + pert;
%  b_noise = randn(length(b_c),1); b_noise = b_noise/norm2(b_noise);
%  b_c = b_c + 1e-12 * norm_b_c * b_noise; 

diagR_c = numtype(zeros(1000,1));
supdR_c = numtype(zeros( 999,1)); % super diagonal

beta_u_c = b_c;
beta_c = norm2(beta_u_c);
u_c = beta_u_c / beta_c;

alpha_v_c = OP'*u_c;
alpha_c = norm2(alpha_v_c);
v_c = alpha_v_c / alpha_c;

diagB_c(i) = alpha_c;

w_c = v_c;
x_c = zeros(n,1);
phi_bar_c = beta_c;
rho_bar_c = alpha_c;

if params.verify
  V(:,i)=v_c;
  U(:,i)=u_c;
  B(i,i)=alpha_c;
  beta_1 = beta_c;
end

% set up the inconsistent LSQR 

b_i = randn(m,1); % possibly inconsistent
norm_b_i = norm2(b_i);

diagR_i = numtype(zeros(1000,1));
supdR_i = numtype(zeros( 999,1)); % super diagonal

beta_u_i = b_i;
beta_i = norm2(beta_u_i);
u_i = beta_u_i / beta_i;

alpha_v_i = OP'*u_i;
alpha_i = norm2(alpha_v_i);
v_i = alpha_v_i / alpha_i;

diagB_i(i) = alpha_i;

w_i = v_i;
x_i = zeros(n,1);
phi_bar_i = beta_i;
rho_bar_i = alpha_i;

while true % LSQR
  % continue the bidiagonalization
  if mod(i,10000)==0; fprintf('starting iteration %d\n',i); end

  %=== consistent iteration ===
  
  beta_u_c = OP*v_c - alpha_c*u_c;
  beta_c = norm2(beta_u_c);
  u_c = beta_u_c / beta_c;
  
  alpha_v_c = OP'*u_c - beta_c*v_c;
  alpha_c = norm2(alpha_v_c);
  v_c = alpha_v_c / alpha_c;
   
  % construct and apply next orthogonal transformation
  rho_c = sqrt( rho_bar_c*rho_bar_c + beta_c*beta_c);
  c_c       = rho_bar_c / rho_c;
  s_c       = beta_c / rho_c;
  theta_c   = s_c*alpha_c;
  rho_bar_c = -c_c*alpha_c;
  phi_c     = c_c*phi_bar_c;
  phi_bar_c = s_c*phi_bar_c;
   
  % update x, w
  x_c = x_c + (phi_c/rho_c)*w_c;
  w_c = v_c - (theta_c/rho_c)*w_c;
  
  % update R, B
  diagR_c(i) = rho_c;
  supdR_c(i) = theta_c;
  
  diagB_c(i+1) = alpha_c;
  subdB_c(i)   = beta_c;
  
  %=== statistics for the consistent problem ===
  if params.trace && mod(i,trace_interval)==0
    if params.trace_with_svd
      sigma_min_R_c(i) = min(svd(diag(diagR_c(1:i))+diag(supdR_c(1:i-1),1)));
    else
      [smaxinvR_c,dummy] = power_bidiag_inverse(diagR_c(1:i),supdR_c(1:i-1));
      sigma_min_R_c(i) = 1/smaxinvR_c;
    end            
  end

  
  %est_v_c(i) = norm(A*v_c) / norm(v_c); 
  
  diff = x_c - x_orig;
  normdiff = norm2(diff);
  normdiff_c(i) = normdiff;
  
  if params.projtrace
      PT(i,1) = abs(diff'*v_pt_min)/normdiff;
      PT(i,2) = abs(diff'*v_pt_mid)/normdiff;
      PT(i,3) = abs(diff'*v_pt_max)/normdiff;
      
      projections(i,:) = abs(VA'*diff);
  end
  
  % normdiff_c(i) = abs( v_min_exact' * (diff/normdiff) );
  % ??? normdiff_c(i) = abs( v_min_exact' * (diff) );
  
  Adiff = A*diff;
  normAdiff = norm2(Adiff);

  normx_c(i) = norm2(x_c);
  % on some matrices with condition number exactly 1, we get diff=0
  % (exactly), normally in the first iteration. This causes the 
  % algorithm to halt, but sigma_min_est=NaN so we never set it.
  % In this case, the initial value of sigma_min=sigma_max
  % is correct and we get a correct result.
  % On any other matrix, the probability of getting diff=0 is zero.
  sigma_min_est  = normAdiff / normdiff;
  sigma_min_d(i) = sigma_min_est;
  normr_est_c(i) = phi_bar_c    / (sigma_max*normx_c(i) + norm_b_c);
  normr_c(i)     = norm2(Adiff) / (sigma_max*normx_c(i) + norm_b_c);
  %ratio = sigma_min_est / normr_c(i);
  if sigma_min_est < sigma_min
   sigma_min = sigma_min_est;
    v_min     = diff;
    %ratio_at_min = ratio;
    iter_at_min = i;
  end
  if sigma_min_est > sigma_min
      estimate_rising = true;
  end
  
  %=== inconsistent iteration ===
  
  beta_u_i = OP*v_i - alpha_i*u_i;
  beta_i = norm2(beta_u_i);
  u_i = beta_u_i / beta_i;
  
  alpha_v_i = OP'*u_i - beta_i*v_i;
  alpha_i = norm2(alpha_v_i);
  v_i = alpha_v_i / alpha_i;
   
  % construct and apply next orthogonal transformation
  rho_i = sqrt( rho_bar_i*rho_bar_i + beta_i*beta_i);
  c_i       = rho_bar_i / rho_i;
  s_i       = beta_i / rho_i;
  theta_i   = s_i*alpha_i;
  rho_bar_i = -c_i*alpha_i;
  phi_i     = c_i*phi_bar_i;
  phi_bar_i = s_i*phi_bar_i;
   
  % update x, w
  x_i = x_i + (phi_i/rho_i)*w_i;
  w_i = v_i - (theta_i/rho_i)*w_i;
  
  % update R, B
  diagR_i(i) = rho_i;
  supdR_i(i) = theta_i;
  
  diagB_i(i+1) = alpha_i;
  subdB_i(i)   = beta_i;

  %=== statistics for the inconsistent problem ===
  if params.trace && mod(i,trace_interval)==0
    if params.trace_with_svd
      sigma_min_R_i(i) = min(svd(diag(diagR_i(1:i))+diag(supdR_i(1:i-1),1)));
    else
      [smaxinvR_i,dummy] = power_bidiag_inverse(diagR_i(1:i),supdR_i(1:i-1));
      sigma_min_R_i(i) = 1/smaxinvR_i;
    end            
  end

  r_i = A*x_i - b_i;
  normx_i(i) = norm2(x_i);
  normr_est_i(i) = phi_bar_i    / (sigma_max*normx_i(i) + norm_b_i);
  normr_i(i)     = norm2(r_i)   / (sigma_max*normx_i(i) + norm_b_i);
  normAtr_i(i)   = norm2(A'*r_i) / (sigma_max*norm2(r_i));
  
  sigma_min_est = norm2(b_i+r_i)/normx_i(i);
  % [sigma_min_est normx_i(i) norm(r_i)]
  sigma_min_x(i) = sigma_min_est;
  
  % params.verify and tracing

  if params.verify
    if i==1; theta_1 = theta_c; end
    V(:,i+1)=v_c;
    U(:,i+1)=u_c;
    if explicit_r
      P(:,i) = p/norm2(p);
    end
  end
    
  i = i+1;

  %=== termination conditions ===
  %class(~params.fixedits)
  %class(normr_c < params.restol)
  %[normr_c(i-1) params.restol]

  % if the matrix is ill conditioned, use a stricter 
  % residual convergence tolerance
  if sigma_max / sigma_min > params.highcond_thresh; 
    params.restol = params.highcond_restol;
  end
  
  if normr_c(i-1) < params.restol; 
      if ~stop
          stop = true;
          stop_res = i-1;
      end
  end
  %if shift==0 && ratio > params.ratiotol && estimate_rising; 
  %    if ~stop
  %        stop = true;
  %        stop_ratio = i-1;
  %    end
  %end
  
  if shift==0 && normdiff_c(i-1) < normdifftol
      if ~stop
          stop = true;
          stop_normdiff = i-1;
      end
  end
  
  if sigma_max / sigma_min > params.cond_thresh; 
      if ~stop
          stop = true;
          stop_cond = i-1;
      end
  end
  
  if ~params.fixedits && ~params.extraits && stop; 
      disp('stopping');
      break; 
  end
  if stop && ~params.fixedits && params.extraits ... 
      && ((~isnan(stop_normdiff) && i-1 > stop_normdiff*1.25) ...
          || (~isnan(stop_res) && i-1 > stop_res*1.25) ... 
          || (~isnan(stop_cond) && i-1 > stop_cond*1.25)); 
      disp('extra its stop');
      break; 
  end
  if (i>=params.maxits); 
      break; 
  end
end

k = i-1; % last step in the algorithm
k 
sigma_min = (sigma_min.^2 - shift^2).^0.5;

% Lanczos estimate
[smaxinvR,dummy] = power_bidiag_inverse(diagR_c(1:k),supdR_c(1:k-1));
sigma_min_R_c_final = 1/smaxinvR
%disp('constructing R');
%R_c = diag(diagR_c(1:k))+diag(supdR_c(1:k-1),1);
%disp('svd(R)');
%sR_c = svd(R_c);
%disp('min(svd(R))');
%sigma_R_c_corrected = (sR_c.^2 - shift^2).^0.5;
%sigma_min_R_c_final = min(double(sigma_R_c_corrected))

fprintf(1,'diff estimate=%.10e at \at iteration=%d, Lanczos estimate=%.2e\n',double(sigma_min), ...
                                                                                      double(iter_at_min), ...
                                                                                      double(sigma_min_R_c_final));

[normr_i(k) normAtr_i(k) normx_i(k)]

min(sigma_min_x(1:k))
    
stats.iterations = k;
stats.sigma_min_R = sigma_min_R_c_final;
if isnan(stop_normdiff) && isnan(stop_res) && isnan(stop_cond)
  stats.converged = false;
else
  stats.converged = true;
end
stats.stop_error    = stop_normdiff;
stats.stop_residual = stop_res;
stats.stop_cond     = stop_cond;

if params.verify
  R_c = diag(diagR_c(1:k))+diag(supdR_c(1:k-1),1);
  B_c = diag(diagB_c(1:k))+diag(subdB_c(1:k-1),-1);
  B_c(k+1,k) = subdB_c(k);
  
  R_i = diag(diagR_i(1:k))+diag(supdR_i(1:k-1),1);
  B_i = diag(diagB_i(1:k))+diag(subdB_i(1:k-1),-1);
  B_i(k+1,k) = subdB_i(k);

  %size(B)
  
  fprintf(1,'Performed %d iterations\n',k);

  fprintf(1,'Relative Residual=%.0e\n',norm2(A*x_c-b_c)/(sA(1)*norm2(b_c)));

  fprintf(1,'Residual in (3.2)=%.0e\n',norm2(U(:,1)*beta_1 - b_c));
  fprintf(1,'Residual in (3.3)=%.0e\n',norm2(A*V(:,1:k) - U(:,1:k+1)*B_c(1:k+1,1:k)));
  fprintf(1,'Residual in (3.4)=%.0e\n',norm2(A'*U(:,1:k+1) - V(:,1:k)*(B_c(1:k+1,1:k))'-alpha_c*v_c*[zeros(1,k) 1]));

  
  %V(:,1:k)'*V(:,1:k)
  fprintf(1,'Orthogonality Residual for V=%.0e\n',norm2(V(:,1:k)'*V(:,1:k) - eye(k)));
  fprintf(1,'Orthogonality Residual for U=%.0e\n',norm2(U(:,1:k+1)'*U(:,1:k+1) - eye(k+1)));

  if explicit_r
    fprintf(1,'Orthogonality Residual for P=%.0e\n',norm2(P(:,1:k)'*P(:,1:k) - eye(k)));
    fprintf(1,'Residual in (3.6)=%.0e\n',norm2(V(:,1)*theta_1 - A'*b));
    fprintf(1,'Residual in (3.7)=%.0e\n',norm2(A*V(:,1:k) - P(:,1:k)*R(1:k,1:k)));
  end
  
  sR_c = svd(R_c(1:k,1:k));
  sB_c = svd(B_c(1:k+1,1:k));

  sR_i = svd(R_i(1:k,1:k));
  sB_i = svd(B_i(1:k+1,1:k));

  fprintf(1,'Bound (5.9): %.3e == %0.3e <= %.3e\n',1/sR_c(k),1/sB_c(k),1/sA(n));

end
  
%if params.verify
%  sigma_min_R = (sigma_min_R.^2 - shift^2).^0.5;
%  sigma_min_B = (sigma_min_B.^2 - shift^2).^0.5;
%  sigma_min_d = (sigma_min_d.^2 - shift^2).^0.5;
%   
%  fprintf(1,'sigma_min(A)=%.2e, diff estimate=%.2e, relative error=%.0e\n',sA(n),sigma_min,abs(sA(n)-sigma_min)/sA(n));
%end

close all;
if params.trace
    
    k = i-1;
    if false
    subplot(2,2,1);
    semilogy(1:k,normr_c(1:k),'bs',1:k,normr_est_c(1:k),'+',1:k,normx_c(1:k),'bx');
    maxerr = max(abs( (norm_r(1:k) - norm_r_est(1:k))./norm_r(1:k)));
    title(sprintf('consistent: max relative error between phibar and norm(r) = %.0e',maxerr));
    legend('norm(r)','phibar \approx norm(r)','1/norm(x)');
    xlabel('iteration');
    a = axis;
    data = [normr_c(1:k) normr_est_c(1:k) normx_c(1:k)];
    axis([ a(1) a(2) min(eps,min(min(data)))/3 max(max(data))*3]);
    h=line( [a(1) a(2)] ,[ eps eps ] );
    set(h,'Color',[1 0 0]); % make the machine espilon line red
    end
    
    %subplot(2,2,1);
    figure(1);
    hh1 = semilogy(1:k,normr_c(1:k),'-bs',1:k,sigma_min_R_c(1:k),'-ro', 1:k,sigma_min_d(1:k),'-m>',1:k,normdiff_c(1:k),'-k^');  
    hold on;
    hl = legend('$\Vert r^{(t)} \Vert / \left( \hat{\sigma}_{\max} \Vert x^{(t)} \Vert + \Vert b \Vert \right)$', ...
                '$\sigma_{\min}(R^{(t)})$', ...
                '$\Vert A d^{(t)} \Vert / \Vert d^{(t)} \Vert$','$\Vert d^{(t)} \Vert$','Location','SouthWest');
    set(hl, 'Interpreter', 'latex');     
    set(hh1, 'Visible', 'off');
    hh2 = semilogy(1:k,normr_c(1:k),'-b',1:k,sigma_min_R_c(1:k),'-r', 1:k,sigma_min_d(1:k),'-m',1:k,normdiff_c(1:k),'-k');      
    sel = [1 (1:floor(k/50))*50 k];
    hh3 = semilogy(sel,normr_c(sel),'bs',sel,sigma_min_R_c(sel),'ro', sel,sigma_min_d(sel),'m>',sel,normdiff_c(sel),'k^'); 
    
    xlabel('Iteration');
    ht = title('Random $x^{\star}$ and consistent $b=Ax^{\star}$');
    set(ht, 'Interpreter', 'latex');
    a = axis;
    data = [normr_c(1:k) normr_est_c(1:k) normx_c(1:k)];
    axis([ a(1) a(2) eps/10 10]);
    if ~isnan(stop_res)
        hline = line([stop_res stop_res],[eps/10 10],'Color',[0 1 0]);
        htext = text(stop_res + 0.01*(a(2)-a(1)),0.01,'$\Vert r^{(t)} \Vert$  termination');
        set(htext, 'Interpreter', 'latex');
        set(hline, 'LineWidth', linewidth);
        set(htext, 'FontSize', legend_fontsize);
        
        if (params.extraits)
            stopf = ceil(stop_res * 1.25);   
            hline1 = line([stopf stopf],[eps/10 10],'Color', 'c');
            htext1 = text(stopf + 0.01*(a(2)-a(1)),0.001,'Return');
            set(htext1, 'Interpreter', 'latex');
            set(hline1, 'LineWidth', linewidth);
            set(htext1, 'FontSize', legend_fontsize);
        end
    end
    if ~isnan(stop_normdiff)
        hline = line([stop_normdiff stop_normdiff],[eps/10 10],'Color',[0 1 0]);
        htext = text(stop_normdiff + 0.01*(a(2)-a(1)),0.01,'$\Vert d^{(t)} \Vert$ termination');
        set(htext, 'Interpreter', 'latex');
        set(hline, 'LineWidth', linewidth);
        set(htext, 'FontSize', legend_fontsize);
        
        if (params.extraits)
            stopf = ceil(stop_normdiff * 1.25);   
            hline1 = line([stopf stopf],[eps/10 10],'Color', 'c');
            htext1 = text(stopf + 0.01*(a(2)-a(1)),0.001,'Return');
            set(htext1, 'Interpreter', 'latex');
            set(hline1, 'LineWidth', linewidth);
            set(htext1, 'FontSize', legend_fontsize);
        end
    end
    if ~isnan(stop_cond)
        hline = line([stop_cond stop_cond],[eps/10 10],'Color',[0 1 0]);
        htext = text(stop_cond + 0.01*(a(2)-a(1)),0.01,'$\kappa$ termination');
        set(htext, 'Interpreter', 'latex');
        set(hline, 'LineWidth', linewidth);
        set(htext, 'FontSize', legend_fontsize);
        
        if (params.extraits)
            stopf = ceil(stop_cond * 1.25);   
            hline1 = line([stopf stopf],[eps/10 10],'Color', 'c');
            htext1 = text(stopf + 0.01*(a(2)-a(1)),0.001,'Return');
            set(htext1, 'Interpreter', 'latex');
            set(hline1, 'LineWidth', linewidth);
            set(htext1, 'FontSize', legend_fontsize);
        end
    end
    grid on
    grid minor
    set(gca, 'XMinorGrid', 'off');
    set(gca, 'XGrid', 'on');
   
    set(gca, 'FontSize', fontsize);
    set(ht, 'FontSize', fontsize);
    set(hl, 'FontSize', legend_fontsize);
    set(get(gca, 'XLabel'), 'FontSize', fontsize);
    set(get(gca, 'YLabel'), 'FontSize', fontsize);

    hplots = get(gca, 'Children');
    set(hh1, 'LineWidth', linewidth);
    set(hh1, 'MarkerSize', markersize);
    set(hh2, 'LineWidth', linewidth);
    set(hh3, 'MarkerSize', markersize);

    figure(2);   
    hh1 = semilogy(1:k,normr_i(1:k),'-bs',1:k,normAtr_i(1:k),'-gd',1:k,sigma_min_R_i(1:k),'-ro', 1:k,sigma_min_x(1:k),'-m>');  
    hold on;
    hl = legend('$\Vert r^{(t)} \Vert / \left( \hat{\sigma}_{\max} \Vert x^{(t)} \Vert + \Vert b \Vert \right)$', ...
                '$\Vert A^{T} r^{(t)} \Vert / \left( \hat{\sigma}_{\max} \Vert r^{(t)} \Vert \right)$', ...
                '$\sigma_{\min}(R^{(t)})$', ...
                '$\Vert b + r^{(t)} \Vert / \Vert x^{(t)} \Vert$', 'Location','SouthWest');
    set(hl, 'Interpreter', 'latex');
    set(hh1, 'Visible', 'off');
    hh2 = semilogy(1:k,normr_i(1:k),'-b',1:k,normAtr_i(1:k),'-g',1:k,sigma_min_R_i(1:k),'-r', 1:k,sigma_min_x(1:k),'-m');  
    sel = [1 (1:floor(k/50))*50 k];
    hh3 = semilogy(sel,normr_i(sel),'bs',sel,normAtr_i(sel),'gd',sel,sigma_min_R_i(sel),'ro', sel,sigma_min_x(sel),'m>');  
    xlabel('Iteration');
    ht = title('Random $b$ (inconsistent system)');
    set(ht, 'Interpreter', 'latex');
    a = axis;
    data = [normr_c(1:k) normr_est_c(1:k) normx_c(1:k)];
    axis([ a(1) a(2) eps/10 10]);
    
    grid on
    grid minor
    set(gca, 'XMinorGrid', 'off');
    set(gca, 'XGrid', 'on');
   
    set(gca, 'FontSize', fontsize);
    set(ht, 'FontSize', fontsize);
    set(hl, 'FontSize', legend_fontsize);
    set(get(gca, 'XLabel'), 'FontSize', fontsize);
    set(get(gca, 'YLabel'), 'FontSize', fontsize);

    hplots = get(gca, 'Children');
    set(hh1, 'LineWidth', linewidth);
    set(hh1, 'MarkerSize', markersize);
    set(hh2, 'LineWidth', linewidth);
    set(hh3, 'MarkerSize', markersize);

    %subplot(2,2,3);
    %semilogy(1:k,normr_i(1:k),'bs',1:k,normr_est_i(1:k),'+',1:k,normx_i(1:k),'bx');
    %maxerr = max(abs( (norm_r(1:k) - norm_r_est(1:k))./norm_r(1:k)));
    %title(sprintf('inconsistent: max relative error between phibar and norm(r) = %.0e',maxerr));
    %legend('norm(r)','phibar \approx norm(r)','norm(x)');
    %xlabel('iteration');
    %a = axis;
    %data = [normr_i(1:k) normr_est_i(1:k) normx_i(1:k)];
    %axis([ a(1) a(2) min(eps,min(min(data)))/3 max(max(data))*3]);
    %h=line( [a(1) a(2)] ,[ eps eps ] );
    %set(h,'Color',[1 0 0]); % make the machine espilon line red

    %subplot(2,2,3);
    %semilogy(1:k,sigma_min_d(1:k),'go',1:k,normr_c(1:k),'bs');
    %title(sprintf('\\sigma_{min}(A)=%.0e, min(\\sigma_{min}(B))=%.0e, min(ours)=%.0e', ...
    %              sA(n),min(sigma_min_B(1:k)),min(sigma_min_d(1:k))));
    %legend('our estimate','norm(r)');
    %xlabel('iteration');  
    
    %subplot(2,2,[3 4]);
    if (params.verify)
        figure(3);
        %sA = svd(A);
        %sO = svd(A_orig);
        %sR = svd(R);
        %sR = (sR.^2 - shift^2).^0.5;
        %sA = (sA.^2 - shift^2).^0.5;
        %[min(sO) max(sO) min(sA) max(sA) min(sR) max(sR)]
        h=semilogx(sA,1*ones(length(sA),1),'k+',sR_c,2*ones(length(sR_c),1),'ro',min(sigma_min_d(1:k)),3,'mo',sR_i,4*ones(length(sR_i),1),'rs',min(sigma_min_x(1:k)),5,'ms');
        set(gca,'YTick',[]);
        set(gca,'YTickLabel',[]);
        a = axis;
        axis([ a(1)/2 a(2)*2 0 5]);
        min(sigma_min_d(1:k))
        title(sprintf('exact and estimated singular values'));
        legend('\Sigma(A)','\Sigma(R) (random x)','||A(x-x_*|| / ||x-x_*|| (random x)','\Sigma(R) (random b)','||b + r|| / ||x|| (random b)','Location','NorthWestOutside');
        xlabel('spectrum');  
        a = axis;
        axis([ a(1) a(2) 0 6]);
    end
    
    if params.projtrace
        figure;
        semilogy(1:size(PT,1),PT(:,1),'x',1:size(PT,1),PT(:,2),'o',1:size(PT,1),PT(:,3),'s');
        legend('v min','v middle','v max');
    end
   
else
    figure(1);
    semilogy(1:k,normr_c(1:k),'bs',1:k,sigma_min_d(1:k),'mo',1:k,normdiff_c(1:k),'k^');  
    legend('||r|| / ( ||A|| ||x|| + ||b|| )','||A(x-x*)|| / ||x-x*||','||x-x*||','Location','SouthWest');
    xlabel('iteration');
    title('random x_* and consistent b=Ax*');
    a = axis;
    %data = [normr_c(1:k) normr_est_c(1:k) normx_c(1:k)];
    axis([ a(1) a(2) eps/10 10]);
    if ~isnan(stop_res)
        line([stop_res stop_res],[eps/10 10],'Color',[0 1 0]);
        text(stop_res + 0.01*(a(2)-a(1)),1,'residual-norm termination');
    end
    if ~isnan(stop_normdiff)
        line([stop_normdiff stop_normdiff],[eps/10 10],'Color',[0 1 0]);
        text(stop_normdiff + 0.01*(a(2)-a(1)),1,'||x-x*|| termination');
    end


end
% spy(R)

if params.projtrace
    close all;
    figure
    subplot(3,1,[1 2]);
    imagesc(log10(abs(projections)),[-16 0]);
    %colorbar;
    ylabel('Iteration');
    xlabel('Index of singular triplet');
    ht = title('log_{10} of the projection of the error on singular vectors');
    set(get(gca, 'XLabel'), 'FontSize', fontsize);
    set(get(gca, 'YLabel'), 'FontSize', fontsize);
    set(ht, 'FontSize', fontsize);
    
    subplot(3,1,3);
    semilogy(1:n,diag(SA),'x');
    xlabel('Index of singular triplet');
    ylabel('Singular value');
    set(get(gca, 'XLabel'), 'FontSize', fontsize);
    set(get(gca, 'YLabel'), 'FontSize', fontsize);
    a=axis;
    axis([a(1) a(2) min(diag(SA))*0.1 max(diag(SA))*10]);
end

    function [ssigma_max,vv_max] = power_rect_explicit(A)
      [mm nn] = size(A);  
      vv = numtype(randn(nn,1));
      vv = vv / norm2(vv);
      itcount = klein_lu_power_bound(min(size(A)),0.1,1e-12);
      for k=1:itcount
        vv_max = vv;
        uu = A*vv;
        vv = A'*uu;
        ssigma_max = norm2(vv);
        vv = vv / ssigma_max;
      end
      vv_max = vv;
      ssigma_max = norm2(A*vv) / norm2(vv);
    end

    function [ssigma_max,vv_max] = power_bidiag_inverse(main,sup) % main and super diagonal
      nn = length(main);
      %%RR = diag(main)+diag(sup,1);
      %if size(main,1)>1; main=main'; end
      %if size(sup ,1)>1; sup =sup' ; end
      vv = numtype(randn(nn,1));
      vv = vv / norm2(vv);
      %I = [ 1:n 1:n-1 ];
      %J = [ 1:n 2:n   ];
      %V = [ main sup ];
      %A=sparse(I,J,V,n,n);
      %opts.UT = true;
      %size(main)
      %size(sup)
      %AA = diag(main)+diag(sup,1);
      itcount = klein_lu_power_bound(min(size(A)),0.1,1e-12);
      for kk=1:itcount
        % now do a solve with A and A', which multiplies by the inverse
        % back solve with A
        vv_start = vv;
        
        uu = vv; 
        for ii=nn:-1:1
          uu(ii) = uu(ii) / main(ii);
          if ii>1; uu(ii-1) = uu(ii-1) - sup(ii-1)*uu(ii); end;
        end

        %%uu = RR\vv; % redo it
        
        % forward solve with A'
        vv = uu; 
        for ii=1:nn
          vv(ii) = vv(ii) / main(ii);
          if ii<nn; vv(ii+1) = vv(ii+1) - sup(ii)*vv(ii); end;
        end
        
        %%vv = RR'\uu; % redo it
        
        norm_vv = norm2(vv);
        vv = vv / norm_vv;
      end
      vv_max = vv;
      uu = vv; 
      for ii=nn:-1:1
        uu(ii) = uu(ii) / main(ii);
        if ii>1; uu(ii-1) = uu(ii-1) - sup(ii-1)*uu(ii); end;
      end
      
      %%uu = RR\vv;
      
      ssigma_max = norm2(uu);
      %result = 1/ssigma_max
      %sRR = svd(RR);
      %minmax_svd_R = [min(sRR) max(sRR)] 
    end

    function nn=norm2(xx)
        if strcmp(class(xx),'double')
            nn = norm(xx);
        else
            %nn = sqrt(xx'*xx);
            nn = norm(double(xx));
        end
    end

    function xx = numtype(zz)
        if strcmp(params.type,'sym')
            xx = vpa(zz);
        else
            xx = zz;
        end
    end

   function kl_itcount = klein_lu_power_bound(kl_n,kl_epsilon,kl_delta)
     kl_itcount = ceil( (1/kl_epsilon) * (log(1*n) + log(1/(kl_epsilon*kl_delta))) );
   end
end
