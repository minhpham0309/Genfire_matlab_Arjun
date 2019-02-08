function indarr = My_find_symm_indarr(asize)

    if mod(asize,2) == 1
        CEN = (asize+1)/2;
        indarr = -1*(CEN-1):1:(CEN-1);
    else
        CEN = asize/2;
        indarr = -1*CEN:1:(CEN-1);
    end
    
end