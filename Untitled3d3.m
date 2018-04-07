addpath(genpath(pwd))
format long
a = 32;
ms = a;
mt = a;
mu = a;
n1 = a;
n2 = a;
n3 = a;
nj = n1*n2*n3;
iflag=-1;
eps=1e-12;
xj=zeros(nj,1);
yj=zeros(nj,1);
zj=zeros(nj,1);
cj=zeros(nj,1);
for k3 = -n3/2:(n3-1)/2
 for k2 = -n2/2:(n2-1)/2
    for k1 = -n1/2:(n1-1)/2
       j =  (k1+n1/2+1) + (k2+n2/2)*n1 + (k3+n3/2)*n1*n2;
       xj(j) = pi*cos(-pi*k1/n1);
       yj(j) = pi*cos(-pi*k2/n2);
       zj(j) = pi*cos(-pi*k3/n3);
       cj(j) = sin(pi*j/n1)+1i*cos(pi*j/n2);
    end
 end
end
%xj=pi*rand(nj,1);
%yj=pi*rand(nj,1);
%cj=pi*rand(nj,1);

nk = ms*mt*mu;
sk=zeros(nk,1);
tk=zeros(nk,1);
uk=zeros(nk,1);
for k1 = 1: nk
   sk(k1) = 10*(cos(k1*pi/nk-pi/2));
   tk(k1) = 8*(sin(k1*pi/nk));
   uk(k1) = 10*(cos(k1*pi/nk-pi/2));
end
x=[xj(:) yj(:) zj(:)];
k=[sk(:) tk(:) uk(:)];
fftconst = iflag*1i/ms*2*pi;
ratiofun = @(k,x)exp(fftconst*(k-round(k))*x');
[U,V] = lowrank(k,x/2/pi*ms,ratiofun,eps,2500,2500);

ksub = mod(round(k),ms)+1;
kksub = sub2ind([ms ms ms],ksub(:,1),ksub(:,2),ksub(:,3));
nufft3fun = nufft3III(k,x/2/pi*ms,iflag,ms,500,eps);
fhat1=nufft3fun(cj);
%fhat1=fhat1/nj;

%xsub = mod(round(x/2/pi*ms),ms)+1;
%xxsub = sub2ind([ms ms],xsub(:,1),xsub(:,2));
%spPerm = sparse(xxsub,1:ms^2,ones(1,ms^2),ms^2,ms^2);
r = size(V,2)
[n,ncol] = size(cj);
M = repmat(conj(V),1,ncol).*reshape(repmat(cj,r,1), n, r*ncol);

    ratiofun = @(k,x)exp(fftconst*k*(x-round(x))');
    [k1,k2,k3] = ndgrid(0:ms-1,0:ms-1,0:ms-1);
    k=[k1(:) k2(:) k3(:)];
    [U1,V1] = lowrank(k,x/2/pi*ms,ratiofun,eps,2500,2500);
    size1=size(U1)
    size2=size(V1)
    xsub = mod(round(x/2/pi*ms),ms)+1;
    xxsub = sub2ind([ms ms ms],xsub(:,1),xsub(:,2),xsub(:,3));
    spPerm = sparse(xxsub,1:nj,ones(1,nj),ms^3,nj);
    r1 = size(V1,2)
    [n,ncol] = size(M);
    M1 = repmat(conj(V1),1,ncol).*reshape(repmat(M,r1,1),n,r1*ncol);
    MM1 = reshape(spPerm*M1,ms,ms,ms,r1*ncol);
    MMM1 = fft3(MM1);
    %MMMM1=fftshift(fftshift(MMM1,1),2);
    MM = squeeze( sum( reshape(repmat(U1,1,ncol).*reshape(MMM1,ms^3,r1*ncol), n, r1, ncol), 2) );
    
    
MMM = MM(kksub,:);

fhat = sum(U.*MMM,2);

%fhat = fhat/nj;

[fk,ier]=nufft3d3(nj,xj,yj,zj,cj,iflag,eps,nk,sk,tk,uk);


fk=fk(:);
%error1 = norm(fhat-fhat1,2)/norm(fhat,2)
error1 = norm(fhat-fk)/norm(fk)
error2 = norm(fhat1-fk)/norm(fk)

fid=fopen('./nufftQY/U3r3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U.');
fclose(fid);
fid=fopen('./nufftQY/V3r3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V.');
fclose(fid);
fid=fopen('./nufftQY/U3i3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U.');
fclose(fid);
fid=fopen('./nufftQY/V3i3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V.');
fclose(fid);
fid=fopen('./nufftQY/UU3r3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U1.');
fclose(fid);
fid=fopen('./nufftQY/VV3r3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V1.');
fclose(fid);
fid=fopen('./nufftQY/UU3i3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U1.');
fclose(fid);
fid=fopen('./nufftQY/VV3i3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V1.');
fclose(fid);
%fid=fopen('cj3r.txt','w');
%fprintf(fid,'%12.16f\r\n',cj);
%fclose(fid);
%fid=fopen('cj3i.txt','w');
%fprintf(fid,'%12.16f\r\n',-1i*cj);
%fclose(fid);
%fid=fopen('Re3r3.txt','w');
%fprintf(fid,'%12.16f\r\n',fhat);
%fclose(fid);
%fid=fopen('Re3i3.txt','w');
%fprintf(fid,'%12.16f\r\n',-1i*fhat);
%fclose(fid);
%fid=fopen('xsub33.txt','w');
%fprintf(fid,'%12.16f\r\n',xsub);
%fclose(fid);
%fid=fopen('MM3r3.txt','w');
%fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',MM);
%fclose(fid);
%fid=fopen('MM3i3.txt','w');
%fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*MM);
%fclose(fid);
%fid=fopen('MMM3r3.txt','w');
%fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',MMM);
%fclose(fid);
%fid=fopen('MMM3i3.txt','w');
%fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*MMM);
%fclose(fid);
%fid=fopen('MMMM3r3.txt','w');
%fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',MMMM);
%fclose(fid);
%fid=fopen('MMMM3i3.txt','w');
%fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*MMMM);
%fclose(fid);
