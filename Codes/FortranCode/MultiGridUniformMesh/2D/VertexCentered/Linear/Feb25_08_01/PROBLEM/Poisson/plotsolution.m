function[] = cntrplot2d(var)

IN = './OUT/solution.dat'

clf reset;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.5 0.5 4 4]);

theend = logical(0);
f = fopen(IN,'r');
hold on;
ipatch=0;

[time,count] = fscanf(f, '%f', 1);
[nrvars,count]  =  fscanf(f, '%d', 1);
[dx(1),count] =  fscanf(f, '%f', 1);
[dx(2),count] =  fscanf(f, '%f', 1);

[xl(1),count] =  fscanf(f, '%f', 1);
[xl(2),count] =  fscanf(f, '%f', 1);

[xu(1),count] =  fscanf(f, '%f', 1);
[xu(2),count] =  fscanf(f, '%f', 1);

[n(1),count] =  fscanf(f, '%d', 1);
[n(2),count] =  fscanf(f, '%d', 1);

xu = xl+dx.*n;

xlg(1) = xl(1);
xug(1) = xu(1);
xlg(2) = xl(2);
xug(2) = xu(2);
   
[A]=fscanf(f,'%f', [(2+nrvars),(n(1)+1)*(n(2)+1)]); % ghost layer included.
   
drawcont(A,n,var,xl,xu,dx,xlg,xug);

text(xug(1)-1.4,xug(2)+0.2,['time = ' num2str(time)], 'FontSize', 14);

%grid on
t = (xug-xlg)./4;
set(gca,'XTick',xlg(1):t(1):xug(1))
set(gca,'YTick',xlg(2):t(2):xug(2))

ipatch
fclose(f);

function drawcont(A,n,var,xl,xu,dx,xlg,xug)
for j=1:n(2)+1
for i=1:n(1)+1
  x(i,j) = A(    1,(j-1)*(n(1)+1)+i);
  y(i,j) = A(    2,(j-1)*(n(1)+1)+i);
  u(i,j) = A(2+var,(j-1)*(n(1)+1)+i);
end;
end;

colormap('default')
surf(x,y,u,'LineStyle','none')
shading interp;
axis([xlg(1),xug(1),xlg(2),xug(2)])
axis([xlg(1),xug(1),xlg(2),xug(2)])
colorbar




