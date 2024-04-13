//This function calculates boundary indexes of clusters of 1’s.
//Calling Sequence
//c = clustersegment(s)
//Parameters 
//s: scalar, vector or matrix of real numbers (clusters of 1s) 
//c: output variable, cell array of size 1 by N, where N is the number of rows in s
//Description
//This function calculates boundary indexes of clusters of 1’s.
//This function calculates the initial and end indices of the sequences of 1's present in the input argument.
//The output variable c is a cell array of size 1 by N, where N is the number of rows in s and each element has two rows indicating the initial index and end index of the cluster of 1's respectively. The indexing starts from 1.
//Examples
//y = clustersegment ([0,1,0,0,1,1])
//y  =
//    2.    5.  
//    2.    6.  


funcprot(0)
function contRange = clustersegment(xhi)
    bool_discon = diff(xhi, 1, 2);
    [Np, Na] = size(xhi);
    contRange = cell(1, Np);

    for i = 1:Np
        idxUp = find(bool_discon(i,:) > 0) + 1;
        idxDwn = find(bool_discon(i,:) < 0);
        tLen = length(idxUp) + length(idxDwn);

        if xhi(i,1) == 1
            contRange{i}(1) = 1;
            contRange{i}(2:2:tLen+1) = idxDwn;
            contRange{i}(3:2:tLen+1) = idxUp;
        else
            contRange{i}(1:2:tLen) = idxUp;
            contRange{i}(2:2:tLen) = idxDwn;
        end

        if xhi(i,Na) == 1
            contRange{i}(length(contRange{i})+1) = Na;
        end

        tLen = length(contRange{i});
        if tLen ~= 0
            new_contRange = zeros(2, tLen / 2);
            for j = 1:tLen / 2
                if ~isempty(contRange{i})
                    new_contRange(1, j) = contRange{i}(2*j-1);
                    new_contRange(2, j) = contRange{i}(2*j);
                end
            end
            contRange{i} = new_contRange;
        end
        disp(contRange{i});
    end

    if Np == 1
        contRange = cell2mat(contRange);

    end
end

//Test Cases
//xhi1 = [0 0 1 1 1 0 0 1 0 0 0 1 1];
//xhi2 = [1,1,1,1,1,1,1,1,1,1,1,1,1];
//xhi3 = [0,0,0,0,0,0,0,0,0,0,0,0,0];
//xhi4 = [
//  0 0 1 1 1 0 0 1 0 0 0 1 1;
//  1 0 1 1 1 0 0 1 0 0 0 1 1;
//  1 1 0 0 1 1 0 0 0 1 0 0 1
//];
//
//Call the function and assign the output to variables
//y1 = clustersegment(xhi1);
//y2 = clustersegment(xhi2);
//y3 = clustersegment(xhi3);
//y4 = clustersegment(xhi4);


