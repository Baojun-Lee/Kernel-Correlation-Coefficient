function e=gauss_kernel(x,y,u)
    e=exp(-(x-y).^2/(u*u));
    