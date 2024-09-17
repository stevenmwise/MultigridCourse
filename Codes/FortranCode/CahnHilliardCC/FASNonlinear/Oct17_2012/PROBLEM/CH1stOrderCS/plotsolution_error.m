function[l2error] = cntrplot2d(nn,var,g1,g2)

s1 = ['0000000' num2str(nn)];
s2 = s1((length(s1)-3):length(s1));

dir1 =['./OUT' num2str(g1) '/']

IN1  = [dir1 'm' s2 '.dat']

clf reset;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.5 0.5 4 4]);

theend = logical(0);
f = fopen(IN1,'r');
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

n1 = n;

xu = xl+dx.*n;

xlg(1) = xl(1);
xug(1) = xu(1);
xlg(2) = xl(2);
xug(2) = xu(2);
   
[A1]=fscanf(f,'%f', [(2+nrvars),n(1)*n(2)]); % ghost layer included.

ipatch
fclose(f);

dir2 =['./OUT' num2str(g2) '/']

IN2  = [dir2 'm' s2 '.dat']

clf reset;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.5 0.5 4 4]);

theend = logical(0);
f = fopen(IN2,'r');
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

n2 = n;

xu = xl+dx.*n;

xlg(1) = xl(1);
xug(1) = xu(1);
xlg(2) = xl(2);
xug(2) = xu(2);
   
[A2]=fscanf(f,'%f', [(2+nrvars),n(1)*n(2)]); % ghost layer included.

ipatch
fclose(f);

if n1(1) ~= 2*n2(1)
  'size error'
  return
end
   
[l2error] = drawcont(A1,A2,n,var,xl,xu,dx,xlg,xug);

text(xug(1)-1.4,xug(2)+0.2,['time = ' num2str(time)], 'FontSize', 14);

%grid on
t = (xug-xlg)./4;
set(gca,'XTick',xlg(1):t(1):xug(1))
set(gca,'YTick',xlg(2):t(2):xug(2))

function[l2error] = drawcont(A1,A2,n,var,xl,xu,dx,xlg,xug)
for j=1:n(2)
for i=1:n(1)
  x(i,j) = A2(1, (j-1)*n(1)+i);
  y(i,j) = A2(2, (j-1)*n(1)+i);
  u1(i,j) = 0.25*(A1(2+var,((2*j-1)-1)*2*n(1)+(2*i-1)) ...
          +       A1(2+var,((2*j-1)-1)*2*n(1)+(2*i  )) ...
          +       A1(2+var,((2*j  )-1)*2*n(1)+(2*i-1)) ...
          +       A1(2+var,((2*j  )-1)*2*n(1)+(2*i  )));
  u2(i,j) = A2(2+var,(j-1)*n(1)+i);
  e(i,j) = u1(i,j)-u2(i,j);
  e2(i,j) = e(i,j)*e(i,j);
end;
end;

format long
l2error = sqrt(sum(sum(e2,'double'),'double')/(n(1)*n(2)));

colormap('default')
surf(x,y,e,'LineStyle','none')
shading interp;
axis([xlg(1),xug(1),xlg(2),xug(2)])
axis([xlg(1),xug(1),xlg(2),xug(2)])
colorbar




