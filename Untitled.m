    N   = 3;
    mu    = 0.72;
    rho   = 0.9895;
    sigma = 0.0034;
    
    delta = 0.0081;
    alpha = 0.72;
    A     = 0.158;
    lambda=0.72;
    b     = 0.4;
    kappa = 0.6;
    beta  = 0.999;
    
    tmp0 = 0;
    irow = 1;
    for j = 1 : N
        tmp0 = tmp0 + PP(irow, j) *((1-lambda)*(grid(j)-b) - kappa*lambda*grid(j)  ...
               * exp(xc(j)) + (1-delta)*kappa / qx(j));
    end 