funcprot(0)
function [SOS, G] = zp2sos(z, p, k, DoNotCombineReal)
//This function converts filter poles and zeros to second-order sections.
//Calling Sequence
//[sos] = zp2sos(z)
//[sos] = zp2sos(z, p)
//[sos] = zp2sos(z, p, k)
//[sos, g] = zp2sos(...)
//Parameters 
//z: column vector
//p: column vector
//k: real or complex value, default value is 1
//Description
//This function converts filter poles and zeros to second-order sections.
//The first and second parameters are column vectors containing zeros and poles. The third parameter is the overall filter gain, the default value of which is 1.
//The output is the sos matrix and the overall gain.
//If there is only one output argument, the overall filter gain is applied to the first second-order section in the sos matrix.
//Examples
//zp2sos([1, 2, 3], 2, 6)
//ans =
//    6  -18   12    1   -2    0
//    1   -3    0    1    0    0
    
    if nargin < 3
        k = 1;
    end
    if nargin < 2
        p = [];
    end
    
    if nargin < 4
         DoNotCombineReal = 0;
    end
    
    function [ac, ar] = custom_cplxreal(a)
        tol = 100*(1e-10);
        ac = []; // Initialize empty vector for positive conjugates
        ar = []; // Initialize empty vector for real parts
        n = length(a);
        for i = 1:n-1
            for j = 1:n-i
                if abs(a(j)) > abs(a(j+1))
                    temp = a(j);
                    a(j) = a(j+1);
                    a(j+1) = temp;
                elseif abs(a(j)) == abs(a(j+1)) then
                        if real(a(j)) == real(a(j+1))
                            if imag(a(j)) > imag(a(j+1)) then
                                temp = a(j);
                                a(j) = a(j+1);
                                a(j+1) = temp;
                            end
                        end     
                end
            end
        end
        acp = zeros(1, length(a)); 
        idx_acp = 1;

        for idx = 1:length(a)
            if imag(a(idx)) ~= 0 && idx < length(a) && abs(a(idx) - conj(a(idx+1))) < tol
                acp(idx_acp) = a(idx);
                acp(idx_acp + 1) = conj(a(idx));
                idx_acp = idx_acp + 2;
            else
                acp(idx_acp) = a(idx);
                idx_acp = idx_acp + 1;
            end
        end

        if pmodulo(length(a), 2) == 1
            zcp(length(a)) = a(length(a));
        end
        acp = unique(acp);

        ac1 = acp(find(imag(acp) ~= 0));

        if pmodulo(length(ac1), 2) ~= 0
            error("Some complex numbers could not be paired.");
        end

        for idx = 1:2:length(ac1)-1
            if abs(ac1(idx) - conj(ac1(idx+1))) > tol
                error("Error: All complex numbers are not exact conjugates.");
            end
        end
    
        ac = acp(find(imag(acp) > 0));
        ar = acp(find(imag(acp) == 0));
    end
        
    [zc,zr] = custom_cplxreal(z(:));
    [pc,pr] = custom_cplxreal(p(:));
    
    nzc = length(zc);
    npc = length(pc);

    nzr = length(zr);
    npr = length(pr);
    
    if DoNotCombineReal
        
        SOS = zeros(npc + ceil(npr / 2) + nzc + ceil(nzr / 2), 6);

        for count = 1:npc
            SOS(count, 4:6) = [1, -2 * real(pc(count)), abs(pc(count))^2];
        end

        for count = 1:floor(npr/2)
            SOS(count + npc, 4:6) = [1, -pr(2 * count - 1) - pr(2 * count), pr(2 * count - 1) * pr(2 * count)];
        end

        if modulo(npr, 2) == 1
            SOS(npc + floor(npr / 2) + 1, 4:6) = [0, 1, -pr($)];
        end

        for count = 1:nzc
            SOS(count + npc + ceil(npr / 2), 1:3) = [1, -2 * real(zc(count)), abs(zc(count))^2];
        end

        for count = 1:floor(nzr/2)
            SOS(count + nzc + npc + ceil(npr / 2), 1:3) = [1, -zr(2 * count - 1) - zr(2 * count), zr(2 * count - 1) * zr(2 * count)];
        end

        if modulo(nzr, 2) == 1
            SOS(nzc + floor(nzr / 2) + npc + ceil(npr / 2) + 1, 1:3) = [0, 1, -zr($)];
        end

        if npc + ceil(npr / 2) > nzc + ceil(nzr / 2)
            for count = nzc + ceil(nzr / 2) + 1:npc + ceil(npr / 2)
                SOS(count, 1:3) = [0, 0, 1];
            end
        else
            for count = npc + ceil(npr / 2) + 1:nzc + ceil(nzr / 2)
                SOS(count, 4:6) = [0, 0, 1];
            end
        end
    else
        // Handling complex conjugate poles
        for count = 1:npc
            SOS(count, 4:6) = [1, -2 * real(pc(count)), abs(pc(count))^2];
        end
        //disp(SOS(count, 4:6));

        // Handling pair of real poles
        for count = 1:floor(npr/2)
            SOS(count+npc, 4:6) = [1, - pr(2 * count - 1) - pr(2 * count), pr(2 * count - 1) * pr(2 * count)];
        end
        //disp(SOS(count+npc, 4:6));

        // Handling last real pole (if any)
        if modulo(npr, 2) == 1
            SOS(npc + floor(npr / 2) + 1, 4:6) = [0, 1, -pr($)];
        end

        // Handling complex conjugate zeros
        for count = 1:nzc
            SOS(count, 1:3) = [1, -2 * real(zc(count)), abs(zc(count))^2];
        end
        //disp(SOS(count, 1:3));

        // Handling pair of real zeros
        for count = 1:floor(nzr / 2)
            SOS(count+nzc, 1:3) = [1, - zr(2 * count - 1) - zr(2 * count), zr(2 * count - 1) * zr(2 * count)];
        end
        //disp(SOS(count+nzc, 1:3));

        // Handling last real zero (if any)
        if modulo(nzr, 2) == 1
            SOS(nzc + floor(nzr / 2) + 1, 1:3) = [0, 1, -zr($)];
        end

        // Completing SOS if needed (sections without pole or zero)
        if npc + ceil(npr / 2) > nzc + ceil(nzr / 2)
            for count = nzc + ceil(nzr / 2) + 1 : npc + ceil(npr / 2) // sections without zero
                SOS(count, 1:3) = [0, 0, 1];
            end
        else
            for count = npc + ceil(npr / 2) + 1 : nzc + ceil(nzr / 2) // sections without pole
                SOS(count, 4:6) = [0, 0, 1];
            end
            //disp(SOS(count, 4:6));
        end
    end

    if ~exists('SOS') then
        SOS = [0, 0, 1, 0, 0, 1]; 
    end

    for count = 1:size(SOS, 1)
        B = SOS(count, 1:3);
        A = SOS(count, 4:6);

        while B(1) == 0 && A(1) == 0
            A(1) = [];
            A(3) = 0;
            B(1) = [];
            B(3) = 0;
        end
        SOS(count, :) = [B, A];
    end

    if (nargout < 2)
        SOS(1, 1:3) = k * SOS(1, 1:3);
    else
        G = k;
    end
    
end


//test case
//sos = zp2sos ([]);
//sos = zp2sos ([], []);
//sos = zp2sos ([], [], 2);
//[sos, g] = zp2sos ([], [], 2);
//sos = zp2sos([], [0], 1);
//sos = zp2sos([0], [], 1);
//sos = zp2sos([-1-%i, -1+%i], [-1-2*%i, -1+2*%i], 10);

