dirname = './data';
EMPTY = 999999;
Tday = 86400;
Tref = 0;

fname = [dirname,'/profdata.dat'];

fid = fopen(fname,'rb');
numTotalDataPoints = fread(fid,1,'int32');
numInterpPoints = fread(fid,1,'int32');
Nkmax = fread(fid,1,'int32');
nsteps = fread(fid,1,'int32');
ntoutProfs = fread(fid,1,'int32');  
dt = fread(fid,1,'float64');
dz = fread(fid,Nkmax,'float64');
dataIndices = fread(fid,numTotalDataPoints,'int32');
dataXY = fread(fid,2*numTotalDataPoints,'float64');
xv = reshape(fread(fid,numInterpPoints*numTotalDataPoints,'float64'),numInterpPoints,numTotalDataPoints);
yv = reshape(fread(fid,numInterpPoints*numTotalDataPoints,'float64'),numInterpPoints,numTotalDataPoints);
fclose(fid);

dataX = dataXY(1:2:end);
dataY = dataXY(2:2:end);

dzh = 0.5*(dz(1:end-1)+dz(2:end));
z = zeros(Nkmax,1);
z(1) = -dz(1)/2;
z(2:Nkmax) = z(1) - cumsum(dzh);
dZ = dz*ones(1,length(dataX));

x = sqrt((dataX-dataX(1)).^2+(dataY-dataY(1)).^2)';
X = ones(Nkmax,1)*x;
Z = z*ones(1,numTotalDataPoints);

% There are 5 locations
loc1 = 1;
loc2 = 2;

fname = [dirname,'/u.dat.prof'];
ufid = fopen(fname,'rb');
fname = [dirname,'/T.dat.prof'];
Tfid = fopen(fname,'rb');

component = 1;

udata = fread(ufid,'float64');
nout_u = length(udata)/(3*Nkmax*numInterpPoints*numTotalDataPoints);
udata = reshape(udata,Nkmax,numInterpPoints,numTotalDataPoints,3,nout_u);
u1 = squeeze(udata(:,1,loc1,component,:));
u2 = squeeze(udata(:,1,loc2,component,:));

Tdata = fread(Tfid,'float64');
nout_T = length(Tdata)/(Nkmax*numInterpPoints*numTotalDataPoints);
Tdata = reshape(Tdata,Nkmax,numInterpPoints,numTotalDataPoints,nout_T);
T1 = squeeze(Tdata(:,1,loc1,:));
T2 = squeeze(Tdata(:,1,loc2,:));

T1(find(T1==EMPTY))=nan;
T2(find(T2==EMPTY))=nan;
u1(find(u1==EMPTY))=nan;
u2(find(u2==EMPTY))=nan;

k1 = length(find(~isnan(u1(:,1))));
k2 = length(find(~isnan(u2(:,1))));
d1 = sum(dz(1:k1));
d2 = sum(dz(1:k2));

u10 = sum(u1(1:k1,:).*(dz(1:k1)*ones(1,nout_u)))/d1;
u20 = sum(u2(1:k2,:).*(dz(1:k2)*ones(1,nout_u)))/d2;

t_T = [1:nout_T]*dt*ntoutProfs/Tday;
t_u = [1:nout_u]*dt*ntoutProfs/Tday;

[time_T,Z_T]=meshgrid(t_T,z);
[time_u,Z_u]=meshgrid(t_u,z);

tmax = max(max(time_u));

figure(1);
subplot(4,1,1)
plot(time_u(1,:),u10,'k-');
title('u 1');

subplot(4,1,2)
pcolor(time_u,Z_u,u1)
axis([0 tmax -d1 0]);
shading flat;
colorbar;

subplot(4,1,3)
plot(time_u(1,:),u20,'k-');
title('u 2');

subplot(4,1,4)
pcolor(time_u,Z_u,u2)
axis([0 tmax -d2 0]);
shading flat;
colorbar;

figure(2);
subplot(2,1,1)
contour(time_T,Z_T,T1-Tref*T1(:,1)*ones(1,nout_T),'k-')
axis([0 tmax -d1 0]);
shading flat;
colorbar;
title('T 1');

subplot(2,1,2)
contour(time_T,Z_T,T2-Tref*T2(:,1)*ones(1,nout_T),'k-')
axis([0 tmax -d2 0]);
shading flat;
colorbar;
title('T 2');