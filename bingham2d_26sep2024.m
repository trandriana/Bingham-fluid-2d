clear all
clc
close all

format shortE

global R y0 C eta B r
eta = .25;
r = 3;
rho = 1.75;
B = 2;
R= 1;
C = 16;	
y0 = B/C;
tol = 1.e-5; % Tolerance of the relative error
iterMax = 1000;

M = 50;
N = M;

[u, v, p, div, Du] = stokes2d(eta, r, rho, B, R, C, y0, tol, iterMax, M, N);
x = linspace(0,R, N + 2); h = x(2) - x(1);
xmid = (x(2:end) + x(1:end-1))/2;

plot([0 1 1 0 0], [0 0 1 1 0], 'r')
hold on
 
[X Y] = meshgrid(x(2:end-1), x(2:end-1)');
[Xmid Ymid] = meshgrid(xmid, xmid);
U = (u(1:end-1, :) + u(2:end, :))/2;
V = (v(:, 1:end-1) + v(:, 2:end))/2;
 figure(1)
imagesc(xmid,xmid,Du)

colormap(turbo)
xlim([x(2), x(end-1)])
ylim([x(2), x(end-1)])
% figure('1')
% pcolor(X, Y, region)
 figure(2)
streamslice(Xmid, Ymid, U, V);
% figure('2')
% quiver(Xmid, Ymid, U, V, 'b');

% axis tight
xlim([xmid(1), xmid(end)])
ylim([xmid(1), xmid(end)])


function  [u, v, pre, divu, moduleDu] = stokes2d(eta, r, rho, B, R, C, y0, tol, iterMax, M, N)
    x = linspace(0,R, N+2); h = x(2) - x(1);
     
    xmid = (x(2:end) + x(1:end-1))/2;
     
    A1 = buildA(M+2, N+1);
    A2 = buildA(M+1, N+2);
    
    % Omega1
    u = zeros(M+2, N+1); u0 = u;
    % Omega2
    v = zeros(M+1, N+2); v0 = v;

    % Omega3
    pre = zeros(M, N);
    pre0 = ones(size(pre));
    divu = pre;
    lambda11 = pre; lambda22 = pre;
    
    lambda11 = pre; lambda22 = pre;
    p11 = pre;  p22 = pre;

    % Omega 4
    lambda12 = zeros(M+1, N+1); lambda21 = zeros(M+1, N+1);
    p12 = lambda12; p21 = lambda21;
    
    % Boundary conditions
    ubc = zeros(M, N-1);
    vbc = zeros(M-1, N);
    dL1 = ubc;
    dL2 = vbc;
    
 
%     u(1, :) = 16*(xmid.*(xmid- R)).^2; % boundary coundition BOT
    u(end, :) = 16*(xmid.*(xmid- xmid(end) )).^2; % boundary coundition TOP
%     u(:, 1) = 0.1*ones(size(u(:, 1))); % boundary coundition LEFT
%     u(:, end) = 0.1*ones(size(u(:, end))); % boundary coundition RIGHT
%  
%     v(1, :) = ones(size(v(1, :))); % boundary coundition BOT
%     v(end, :) = ones(size(v(end, :))); % boundary coundition TOP
%     v(:, 1) = 16*(x.*(x- R)).^2; % boundary coundition LEFT
%     v(:, end) = -16*(x.*(x- R)).^2; % boundary coundition RIGHT
    
    ubc(1, :) = r*u(1, 2:end-1)/h^2;
    ubc(end, :) = ubc(end, :) + r*u(end, 2:end-1)/h^2 ;
    ubc(:, 1) = ubc(:, 1) + r*u(2:end-1, 1)/h^2 ;
    ubc(:, end) = ubc(:, end) + r*u(2:end-1, end)/h^2;
    
    vbc(1, :) = r*v(1, 2:end-1)/h^2;
    vbc(end, :) = vbc(end, :) + r*v(end, 2:end-1)/h^2 ;
    vbc(:, 1) = vbc(:, 1) + r*v(2:end-1, 1)/h^2 ;
    vbc(:, end) = vbc(:, end) + r*v(2:end-1, end)/h^2;
    
    F1 = zeros(M, (N-1));
    F2 = zeros((M-1), N);
    
    iter0  = 1;
    res_err0 = 1;
    tic;
    
    dL1 = zeros(M, N-1); dL2 = zeros(M-1, N);
    uu = u;
    uu0 = uu;
    while (res_err0 > tol & iter0 < 1000)
        lambda110 = lambda11; lambda120 = lambda12;
        lambda210 = lambda21; lambda220 = lambda22;
        uu0 = u;
    
        dL1 = diff(lambda11 - r*p11, 1, 2)/h + diff(lambda12(:, 2:end-1) - r*p12(:, 2:end-1), 1, 1)/h;
        dL2 = diff(lambda21(2:end-1, :) - r*p21(2:end-1, :), 1, 2)/h + diff(lambda22 - r*p22, 1, 1)/h;
    %     
        pre = zeros(M, N);
        pre0 = ones(size(pre));
        iter  = 1;
        res_err = 1;
    
        while (res_err > tol & iter < iterMax) 
            pre0 = pre; 
    
            F1 = C - diff(pre, 1, 2)/h + dL1 + ubc;
            F2 = C - diff(pre, 1, 1)/h + dL2 + vbc;
            
            u(2:end-1, 2:end-1) = reshape(-(r*A1/h^2)\F1(:), M, N - 1);
            v(2:end-1, 2:end-1) = reshape(-(r*A2/h^2)\F2(:), M-1, N);
    
            divu = diff(u(2:end-1, :), 1, 2)/h + diff(v(:, 2:end-1), 1, 1)/h;
            pre = pre - rho*divu;
           
            res_err = norm(pre0 - pre)/norm(pre);
            iter = iter + 1;        
        end
        

        D11 = diff(u(2:end-1,:), 1, 2)/h;  D12 = (diff(u, 1, 1) + diff(v, 1, 2))/h/2;
        D21 = (diff(u, 1, 1) + diff(v, 1, 2))/h/2; D22 = diff(v(:, 2:end-1), 1, 1)/h; 
    
        l12  = (lambda12(2:end, 2:end) + lambda12(1:end-1, 2:end)  +  lambda12(2:end, 1:end-1) + lambda12(1:end-1, 1:end-1))/4;
        l21  = (lambda21(2:end, 2:end) + lambda21(1:end-1, 2:end)  +  lambda21(2:end, 1:end-1) + lambda21(1:end-1, 1:end-1))/4;
        DD = (D12(2:end,2:end) + D12(1:end-1,2:end) + D12(2:end,1:end-1) + D12(1:end-1,1:end-1))/4;
    
        X11 =  lambda11 + r*D11; X12 =  l12 + r*DD;
        X21 =  l21 + r*DD;       X22 =  lambda22 + r*D22;
        region = sqrt(X11.^2 + X22.^2 + X12.^2  + X21.^2);
    
        p11 = X11.*(B < region).*(1 - B./(region + (region == 0)))/(eta + r); pp12 = X12.*(B < region).*(1 - B./(region + (region == 0)))/(eta + r);%
        pp21 = X21.*(B < region).*(1 - B./(region + (region == 0)))/(eta + r); p22 = X22.*(B < region).*(1 - B./(region + (region == 0)))/(eta + r);%
    
        p12(2:end-1, 2:end-1) = (pp12(2:end, 2:end) + pp12(1:end-1, 2:end)  +  pp12(2:end, 1:end-1) + pp12(1:end-1, 1:end-1))/4;
        p21(2:end-1, 2:end-1) = (pp21(2:end, 2:end) + pp21(1:end-1, 2:end)  +  pp21(2:end, 1:end-1) + pp21(1:end-1, 1:end-1))/4;
     
        p12(1, :) = D12(1, :); p12(end, :) = D12(end, :); p12(:, 1) = D12(:, 1); p12(:, end) = D12(:, end);
        p21(1, :) = D21(1, :); p21(end, :) = D21(end, :); p21(:, 1) = D21(:, 1); p21(:, end) = D21(:, end);
     
        lambda11 = lambda11  + r*(D11 - p11); lambda12 = lambda12  + r*(D12 - p12);
        lambda21 = lambda21  + r*(D21 - p21); lambda22 = lambda22  + r*(D22 - p22);
    
        res_err0 = norm(lambda110-lambda11)/norm(lambda11) ;
        iter0 = iter0 + 1;
        disp([int2str(iter0),' iterations, ', 'res_err = ', num2str(res_err0), ', elapsed time ', num2str(toc),' sec'])
    end
   
    moduleDu = sqrt(D11.^2 + D22.^2 + 2*DD.^2);
    disp([int2str(iter),' iterations, ', 'res_err = ', num2str(res_err)])

    function A = buildA(M, N)
        I = speye(M - 2, M - 2); 
        II = speye(N - 2, N - 2);
        E = sparse(2:M-2, 1:M-3, 1, M-2, M-2);
        EE = sparse(2:N-2, 1:N-3, 1, N-2, N-2);
        D = E - 4*I + E';
        A = kron(EE, I) + kron(II, D) + kron(EE', I);
    end 


end