function [] = create_plot_2d_data(XC, SOL, ERR, RES, TE)

file = fopen('out/sol.txt','w');
for I=1:length(XC)
  fprintf(file, '%6.12f %6.12f\n', XC(I), SOL(I));
end
fclose(file);

file = fopen('out/err.txt','w');
for I=1:length(XC)
  fprintf(file, '%6.12f %6.12f\n', XC(I), ERR(I));
end
fclose(file);
file = fopen('out/res.txt','w');
for I=1:length(XC)
  fprintf(file, '%6.12f %6.12f\n', XC(I), RES(I));
end
fclose(file);
file = fopen('out/te.txt','w');
for I=1:length(XC)
  fprintf(file, '%6.12f %6.12f\n', XC(I), TE(I));
end
fclose(file);
file = fopen('out/teres.txt','w');
for I=1:length(XC)
  fprintf(file, '%6.12f %6.12f\n', XC(I), RES(I)-TE(I));
end
fclose(file);

end
