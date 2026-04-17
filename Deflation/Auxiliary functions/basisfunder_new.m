function[der] = basisfunder_new(knots, x, mu, p, r)

% basisfunder: Evaluates the r-th derivative of the B-Spline basis functions 
% at x. Based on Algorithm 3.18 of [1].
% INPUT:   
% knots:    Knot vector
% x:        Point(s) in the parametric space (size: N x 1)
% mu:       Knot span number(s) (size: N x 1)
% p:        Spline order
% r:        Derivative order
% OUTPUT:
% der:      Evaluation of the derivative of B-splines at x (size: N x p+1)
%
% Reference:
% [1] T. Lyche and K. Morken. Spline methods draft. University of Oslo, Oslo, 2011.

% Making sure x and mu are both column vectors
if size(x,2)>1
    x=x';
end

if size(mu,2)>1
    mu=mu';
end

N=length(x);
B=basisfun_new(knots, x, mu, p-r);

for k=p-r+1:p

    i1=mu+((1-k):0);
    i2=mu+(1:k);
    t1=reshape(knots(i1(:)), N, k);
    t2=reshape(knots(i2(:)), N, k);
    omega=1./(t2-t1);
    B=[-omega.*B zeros(N,1)]+[zeros(N,1) omega.*B];
end
der=factorial(p)/factorial(p-r)*B;