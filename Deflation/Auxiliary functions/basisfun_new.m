function B = basisfun_new(knots, x, mu, p)

% basisfun:  Evaluates the B-Spline basis functions at x. Based on 
% Algorithm 2.22 of [1].
% INPUT:  
% knots:    Knot vector
% x:        Point(s) in the parametric space (size: N x 1)
% mu:       Knot span number(s) (size: N x 1)
% p:        Spline order
% OUTPUT:
% B:        Evaluation of the B-splines at x (size: N x p+1)
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

Np=length(x);
B=ones(Np,1);

for k=1:p

    i1=mu+((1-k):0);
    i2=mu+(1:k);
    t1=reshape(knots(i1(:)), Np, k);
    t2=reshape(knots(i2(:)), Np, k);
    omega=(x-t1)./(t2-t1);
    B=[(1-omega).*B zeros(Np,1)]+[zeros(Np,1) omega.*B];
end