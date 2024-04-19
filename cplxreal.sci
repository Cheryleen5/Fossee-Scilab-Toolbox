funcprot(0)
function [zc, zr] = cplxreal (z, tol, dim)
//Function to divide vector z into complex and real elements, removing the one of each complex conjugate pair.
//Calling Sequence
//[zc, zr] = cplxreal (z, tol,dim)
//[zc, zr] = cplxreal (z, tol)
//[zc, zr] = cplxreal (z)
//zc = cplxreal (z, tol)
//zc = cplxreal (z)
//Parameters 
//z: vector of complex numbers.
//tol: tolerance for comparisons.
//dim: By default the complex pairs are sorted along the first non-singleton dimension. If dim is specified, then the complex pairs are sorted along this dimension.
//zc: vector containing the elements of z that have positive imaginary parts.
//zr: vector containing the elements of z that are real.
//Description
//Every complex element of z is expected to have a complex-conjugate elsewhere in z. From the pair of complex-conjugates, the one with the negative imaginary part is removed.
//If the magnitude of the imaginary part of an element is less than the tol, it is declared as real.  

  if (nargin < 1 || nargin > 3)
    print_usage ();
  end

  if (isempty (z))
    zc = zeros (size (z));
    zr = zeros (size (z));
    return;
  end

  if typeof(z) == "single" then
    cls = "single";
  else
    cls = "double";
  end

  if (nargin < 2 || isempty (tol))
    tol = 100*(1e-10);
  else
      if (tol >= 1 || tol < 0)
          error("Error: TOL must be a scalar number in the range 0 <= TOL < 1");
      end
  end

  args = cell (1, nargin);
  args{1} = z;
  args{2} = tol;
  if (nargin == 3)
    args{3} = dim;
  else
      dim = 1;
  end

  if nargin == 3
      if dim > ndims(z) && dim ~= 1
        error('Error:Invalid dimension DIM');
      else
        for i = 1:size(z, 1)
            for j = 1:size(z, 2)-1
                for k = 1:size(z, 2)-j
                    if abs(z(i, k)) > abs(z(i, k+1)) then
                        temp = z(i, k);
                        z(i, k) = z(i, k+1);
                        z(i, k+1) = temp;
                    elseif abs(z(i, k)) == abs(z(i, k+1)) then
                        if real(z(i,k)) == real(z(i,k+1))
                            if imag(z(i, k)) > imag(z(i, k+1)) then
                                temp = z(i, k);
                                z(i, k) = z(i, k+1);
                                z(i, k+1) = temp;
                            end
                        end     
                    end
                end
            end
        end
      end

  elseif dim == 1
    n = length(z);
        for i = 1:n-1
            for j = 1:n-i
                if abs(z(j)) > abs(z(j+1))
                    temp = z(j);
                    z(j) = z(j+1);
                    z(j+1) = temp;
                elseif abs(z(j)) == abs(z(j+1)) then
                        if real(z(j)) == real(z(j+1))
                            if imag(z(j)) > imag(z(j+1)) then
                                temp = z(j);
                                z(j) = z(j+1);
                                z(j+1) = temp;
                            end
                        end     
                end
            end
        end
  end
  
    zcp = zeros(1, length(z)); 
    idx_zcp = 1;

    for idx = 1:length(z)
        if imag(z(idx)) ~= 0 && idx < length(z) && abs(z(idx) - conj(z(idx+1))) < tol
            zcp(idx_zcp) = z(idx);
            zcp(idx_zcp + 1) = conj(z(idx));
            idx_zcp = idx_zcp + 2;
        else
            zcp(idx_zcp) = z(idx);
            idx_zcp = idx_zcp + 1;
        end
    end

    if pmodulo(length(z), 2) == 1
        zcp(length(z)) = z(length(z));
    end
    zcp = unique(zcp);

    zc1 = zcp(find(imag(zcp) ~= 0));

    if pmodulo(length(zc1), 2) ~= 0
        error("Some complex numbers could not be paired.");
    end

    for idx = 1:2:length(zc1)-1
        if abs(zc1(idx) - conj(zc1(idx+1))) > tol
            error("Error: All complex numbers are not exact conjugates.");
        end
    end
    
    zc = zcp(find(imag(zcp) > 0));
    zr = zcp(find(imag(zcp) == 0));

endfunction

// Test Case
//[zc, zr] = cplxreal ([]);
//[zc, zr] = cplxreal (1);
//[zc, zr] = cplxreal ([1+%i, 1-%i]);
//[zc, zr] = cplxreal (roots ([1, 0, 0, 1]));
//[zc,zr] = cplxreal (1, 2, 3, 4);
//[zc,zr] = cplxreal (1, ones (2, 3))
//[zc,zr] = cplxreal(1,-1);
//[zc,zr] = cplxreal (1, [], 3)
//[zc,zr] = cplxreal ([[1+2*%i, 1-2*%i, 3-4*%i]  3+4*%i, 7+8*%i, 7-8*%i], [], 2)

