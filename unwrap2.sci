funcprot(0)
function retval = unwrap2(x, tol, dim)
//Unwrap radian phases by adding or subtracting multiples of 2*pi.
//Calling Sequence
//B = unwrap2(X)
//B = unwrap2(X, TOL)
//B = unwrap2(X, TOL, DIM)
//Parameters
//Description
//This function unwraps radian phases by adding or subtracting multiples of 2*pi as appropriate to remove jumps greater than TOL.
//
//    TOL defaults to pi.
//
//Unwrap will work along the dimension DIM.  If DIM is unspecified it defaults to the first non-singleton dimension.
//Examples
//unwrap2([1,2,3])
//ans = 
//        1.    2.    3. 
    if nargin < 1
        error("Not enough input arguments.");
    end

    if ~isreal(x)
        error("X must be real.");
    end

    if nargin < 2 || isempty(tol)
        tol = %pi;
    end

    tol = abs(tol);

    sz = size(x);
    nd = ndims(x);

    if nargin == 3
        if ~(or(type(dim)==[1 5 8]) && isscalar(dim) && dim == fix(dim)) || ~(1 <= dim)
            error("DIM must be an integer and a valid dimension.");
        end
    else
        dim = find(sz > 1, 1);
        if isempty(dim)
            dim = 1;
        end
    end

    rng = 2*%pi;

    if dim > nd || sz(dim) == 1
        retval = x;
        return;
    end

    if ~or(isnan(x(:)) | isinf(x(:)))
        sz(dim) = 1;
        zero_pad = zeros(sz);
        zero_padding = resize_matrix(zero_pad, sz(1),sz(2));
        d = cat(dim, zero_padding, -diff(x, 1, dim));
        p = round(abs(d)./rng) .* rng .* ((d > tol) - (d < -tol));
        retval = cumsum(p, dim) + x;
    else
        if isvector(x)
            retval = x;
            xfin_idx = isinf(x);
            xfin = x(xfin_idx);
            d = cat(dim, 0, -diff(xfin, 1, dim));
            p = round(abs(d)./rng) .* rng .* ((d > tol) - (d < -tol));
            retval(xfin_idx) = xfin + cumsum(p, dim);
        else
            nf_idx = ~isinf(x);
            if all(nf_idx(:))
                retval = x;
                return;
            end
            permuteflag = dim ~= 1;
            if permuteflag
                perm_idx = 1:nd;
                perm_idx([1, dim]) = [dim, 1];
                x = permute(x, perm_idx);
                nf_idx = permute(nf_idx, perm_idx);
                sz([1, dim]) = sz([dim, 1]);
                dim = 1;
            end
            zero_padding = zeros(sz);
            zero_padding(:,dim) = zeros(1, sz(dim));
            x_nf = x(nf_idx);
            x = fill_nonfinite_columnwise(x, nf_idx, zero_padding, sz, nd);
            d = [zero_padding; -diff(x, 1, 1)];
            p = round(abs(d)./rng) .* rng .* ((d > tol) - (d < -tol));
            retval = x + cumsum(p, 1);
            retval(nf_idx) = x_nf;
            if permuteflag
                retval = ipermute(retval, perm_idx);
            end
        end
    end
endfunction

function x = fill_nonfinite_columnwise(x, nonfinite_loc, zero_padding, szx, ndx)
//Replace non-finite values of x, as indicated by logical index nonfinite_loc, with next values.

    flip_idx = 1:ndx;
    flip_idx(1) = szx(1):-1:1;

    nf_front = cumprod(nonfinite_loc, 1);
    nf_back = cumprod(nonfinite_loc(flip_idx{:}), 1)(flip_idx{:});
    nf_middle = nonfinite_loc & ~(nf_back | nf_front);

    locs_before = [diff(nf_middle, 1, 1); zero_padding] == 1;
    locs_after = diff([zero_padding; nf_middle], 1, 1) == -1;
    mid_gap_sizes = find(locs_after) - find(locs_before) - 1;
    x(nf_middle) = repelems(x(locs_after), [1:numel(mid_gap_sizes); mid_gap_sizes'])';

    nf_front = nf_front & ~all(nonfinite_loc, 1);
    front_gap_sizes = (sum(nf_front, 1))(any(nf_front, 1))(:);
    x(nf_front) = repelems(x(locs_after), [1:numel(front_gap_sizes); front_gap_sizes'])';
endfunction

//Test Cases
//x = unwrap2([0;1]);
//x = unwrap2([0;19])
//A = [%pi*(-4), %pi*(-2+1/6), %pi/4, %pi*(2+1/3), %pi*(4+1/2), %pi*(8+2/3), %pi*(16+1), %pi*(32+3/2), %pi*64];
//x = unwrap2(A', %pi, 1);
//x=unwrap2([%inf, 0.5, -1, %nan, %inf, -0.5, 1]); 
//x = unwrap2([%nan,%nan,%nan]);

