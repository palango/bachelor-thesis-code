function [] = create_2d(XC, TC)

file = fopen('out/data.txt','w');
for I=1:length(XC)
  fprintf(file, '%6.12f %6.12f\n', XC(I), TC(I));
end
fclose(file);

end
