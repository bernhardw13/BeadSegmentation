function [Wp,rss] = ellipseproj(W, varargin)
% Projects a set of points onto an ellipse.
% The algorithm is proven to converge and reaches an accuracy of 7-8
% significant digits, and it takes 4-5 iterations per point, on average.
%
% Input arguments:
% W:
%    data points to project to the quadratic curve
%
% Output arguments:
% Wp:
%    data points that lie on the quadratic curve
% rss:
%    residual sum of squares (i.e. the sum of squares of the distances)

% References:
% Nikolai Chernov and H. Ma, "Least squares fitting of quadratic curves
%    and surfaces", Computer Vision, 2011, pages 285--302
% D. Eberly, "Distance from a point to an ellipse in 2D" (2004), Geometric
%    Tools, LLC, http://www.geometrictools.com

% Copyright 2013 Levente Hunyadi

narginchk(2,4);
validateattributes(W, {'numeric'}, {'real','finite','2d'}, mfilename, 'W', 1);
if size(W,1) > size(W,2)
    validateattributes(W, {'numeric'}, {'size',[NaN,2]}, 1);
    transposed = true;
    W = W.';
else
    validateattributes(W, {'numeric'}, {'size',[2,NaN]}, 1);
    transposed = false;
end
if nargin > 2
    narginchk(4,4);
    [center, axes, angle] = varargin{:};
    validateattributes(center, {'numeric'}, {'real','finite','vector','numel',2}, mfilename, 'center', 2);
    validateattributes(axes, {'numeric'}, {'real','finite','vector','numel',2}, mfilename, 'axes', 3);
    validateattributes(angle, {'numeric'}, {'real','finite','scalar'}, mfilename, 'angle', 4);
    center = center(:);
    axes = axes(:);
else
    narginchk(2,2);
    params = varargin{:};
    validateattributes(params, {'numeric'}, {'real','finite','vector','numel',5}, mfilename, 'p', 2);
    params = params(:);
    center = params(1:2);
    axes = params(3:4);
    angle = params(5);
end

% set tolerance
tolerance = 1e-9;

% set maximum number of iterations for vectorized Newton's method
iterations_vectorized = 5;

% set maximum number of iterations for Newton's method
iterations_scalar = 20;

if abs((axes(1)-axes(2)) / axes(1)) < tolerance  % special case: circle
    phi = atan2(W(2,:) - center(2), W(1,:) - center(1));
    Wp = [ axes(1)*cos(phi) ; axes(2)*sin(phi) ];
    Wp = bsxfun(@plus, Wp, center);
else  % ellipse
    a = axes(1);
    b = axes(2);

    % matrix Q for rotating the points and the ellipse to the canonical system
    Q = [cos(angle) sin(angle); -sin(angle) cos(angle)];

    % data points in canonical coordinates
    W = Q * bsxfun(@minus, W, center);
    Wsign = sign(W);
    W = abs(W);

    % set initial value
    t = max( a*W(1,:)-a^2, b*W(2,:)-b^2 );

    % apply Newton' method x(n+1) = x(n)-f(x(n)) / f'(x(n))
    n = 0;  % initialize n (counts iterations)
    while n <= iterations_vectorized
        %f = ( (a*W(1,:))./(t+a^2) ).^2 + ( (b*W(2,:))./(t+b^2) ).^2 - 1;
        %df = -(2*a^2*W(1,:).^2)./(t+a^2).^3 -(2*b^2*W(2,:).^2)./(t+b^2).^3;
        ta2 = t + a^2;
        tb2 = t + b^2;
        pp1 = (a*W(1,:) ./ ta2).^2;
        pp2 = (b*W(2,:) ./ tb2).^2;
        f  = pp1 + pp2 - 1;
        df = -2 * (pp1 ./ ta2 + pp2 ./ tb2);

        % compute next iterate
        u = t - f ./ df;  % x(n+1) = x(n) - f(x(n)) / f'(x(n))

        t = u;
        n = n + 1;
    end

    % apply scalar Newton's method
    tol_a2 = tolerance * a^2;
    for i = 1 : size(W,2)
        n = 0;
        while n <= iterations_scalar
            f = ( (a*W(1,i))/(t(i)+a^2) )^2 + ( (b*W(2,i))/(t(i)+b^2) )^2 - 1;
            if f < 0
                break;
            end

            df = -(2*a^2*W(1,i)^2)/(a^2 + t(i))^3 -(2*b^2*W(2,i)^2)/(b^2 + t(i))^3;
            r = f / df;
            if r < tol_a2
                break;
            end

            % compute next iterate
            u = t(i) - r;  % x(n+1) = x(n) - f(x(n)) / f'(x(n))

            t(i) = u;
            n = n + 1;
        end
    end

    % compute the projection of the point onto the ellipse
    Wp = [(a^2*W(1,:))./(t+a^2) ; (b^2*W(2,:))./(t+b^2)];
    Wp = Wp .* Wsign;

    % rotate back to the original system
    Wp = Q' * Wp;
    Wp = bsxfun(@plus, Wp, center);
end

if nargout > 1
    rss = sum(sum((W - Wp).^2,1));
end

if transposed
    Wp = Wp.';
end
