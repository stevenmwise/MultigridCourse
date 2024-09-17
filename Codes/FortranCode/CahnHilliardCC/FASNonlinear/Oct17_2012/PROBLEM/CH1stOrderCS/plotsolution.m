function[] = cntrplot2d(nn,var)

s1 = ['0000000' num2str(nn)];
s2 = s1((length(s1)-3):length(s1));

dir =['./OUT/']

IN  = [dir 'm' s2 '.dat']
OUT = [dir 'm' s2 '.spinodal.g4.0.jpg']

clf reset;

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
   
[A]=fscanf(f,'%f', [(2+nrvars),n(1)*n(2)]); % ghost layer included.

drawcont(A,n,var,xl,xu,dx,xlg,xug);

t = (xug-xlg)./4;
set(gca,'XTick',xlg(1):t(1):xug(1))
set(gca,'YTick',xlg(2):t(2):xug(2))

axis off

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 4]);
set(gcf, 'PaperPosition', [0.25 0.25 3.5 3.5]);

%print('-djpeg','-r400',OUT)

fclose(f);

function drawcont(A,n,var,xl,xu,dx,xlg,xug)
for j=1:n(2)
for i=1:n(1)
  x(i,j) = A(1, (j-1)*n(1)+i);
  y(i,j) = A(2, (j-1)*n(1)+i);
  u(i,j) = A(2+var, (j-1)*n(1)+i);
end;
end;

colormap('bone')
%caxis([-1 1])
surf(x,y,u,'LineStyle','none')
shading interp;
axis([xlg(1),xug(1),xlg(2),xug(2)])
axis equal
axis([xlg(1),xug(1),xlg(2),xug(2)])
colormap('bone')
%caxis([-1 1])
