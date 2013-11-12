function [] = create_plot_3d_data2(X, Y, Z)

file = fopen('data.txt','w');
for I=1:length(X)
  for J=1:length(Y)
    fprintf(file, '%f %f %f\n', X(I,J), Y(I,J), Z(I,J));
  end
  fprintf(file, '\n');
end
fclose(file);

data = contour(X,Y,Z, 6);%[.1, .3, .5, .7, .9]);
data = data';
save 'contour.txt' data

end
