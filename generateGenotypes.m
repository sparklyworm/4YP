function [genotypes_cell, genotypes] = generateGenotypes(n)
    
    % n = number of genes in genotype
    % eg. n = 4 : 0000, 0001, 0010, .....
    
    genotypes_bin = dec2bin(0:2^n-1);
    genotypes_cell = cellstr(genotypes_bin); % cell array of char
    genotypes = convertCharsToStrings(genotypes_cell); % string array 

end