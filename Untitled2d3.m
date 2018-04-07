addpath(genpath(pwd))
format long
a = 336;
n1=a;
n2=a;
ms=a;
mt=a;
nj=n1*n2;
eps=1e-8;
iflag=-1;
xj=zeros(nj,1);
yj=zeros(nj,1);
cj=zeros(nj,1);
for k1=-n1/2:(n1/2-1)
    for k2=-n2/2:(n2/2-1)
        j=(k2+n2/2+1)+(k1+n1/2)*n2;
        xj(j)=pi*cos(-pi*k1/n1);
        yj(j)=pi*cos(-pi*k2/n2);
        cj(j)=sin(pi*j/n1)+1i*cos(pi*j/n2);
    end
end
nk = ms*mt;
sk=zeros(nk,1);
tk=zeros(nk,1);
for k1 = 1:nk
   sk(k1) = 48*(cos(k1*pi/nk-pi/2));
   tk(k1) = 32*(sin(k1*pi/nk));
end
%xj=pi*rand(nj,1);
%yj=pi*rand(nj,1);
%cj=pi*rand(nj,1);
x=[xj(:) yj(:)];

k=[sk(:) tk(:)];


fftconst = iflag*1i/ms*2*pi;
ratiofun = @(k,x)exp(fftconst*(k-round(k))*x');
[U,V] = lowrank(k,x/2/pi*ms,ratiofun,eps,15000,15000);

ksub = mod(round(k),ms)+1;
kksub = sub2ind([ms ms],ksub(:,1),ksub(:,2));
nufft2fun = nufft2III(k,x/2/pi*ms,iflag,ms,15000,eps);
fhat1=nufft2fun(cj);
%fhat1=fhat1/nj;

%xsub = mod(round(x/2/pi*ms),ms)+1;
%xxsub = sub2ind([ms ms],xsub(:,1),xsub(:,2));
%spPerm = sparse(xxsub,1:ms^2,ones(1,ms^2),ms^2,ms^2);
r = size(V,2)
[n,ncol] = size(cj);
M = repmat(conj(V),1,ncol).*reshape(repmat(cj,r,1), n, r*ncol);

    ratiofun = @(k,x)exp(fftconst*k*(x-round(x))');
    [k1,k2] = ndgrid(0:ms-1);
    k=[k1(:) k2(:)];
    [U1,V1] = lowrank(k,x/2/pi*ms,ratiofun,eps,15000,15000);
    size1=size(U1)
    size2=size(V1)
    xsub = mod(round(x/2/pi*ms),ms)+1;
    xxsub = sub2ind([ms ms],xsub(:,1),xsub(:,2));
    spPerm = sparse(xxsub,1:nj,ones(1,nj),ms^2,nj);
    r1 = size(V1,2)
    [n,ncol] = size(M)
    M1 = repmat(conj(V1),1,ncol).*reshape(repmat(M,r1,1),n,r1*ncol);
    MM1 = reshape(spPerm*M1,ms,ms,r1*ncol);
    MMM1 = fft2(MM1);
    %MMMM1=fftshift(fftshift(MMM1,1),2);
    MM = squeeze( sum( reshape(repmat(U1,1,ncol).*reshape(MMM1,ms^2,r1*ncol), n, r1, ncol), 2) );
    
    
MMM = MM(kksub,:);

fhat = sum(U.*MMM,2);

%fhat = fhat/nj;


fid=fopen('./nufftQY/U2r3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U.');
fclose(fid);
fid=fopen('./nufftQY/V2r3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V.');
fclose(fid);
fid=fopen('./nufftQY/U2i3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U.');
fclose(fid);
fid=fopen('./nufftQY/V2i3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V.');
fclose(fid);
fid=fopen('./nufftQY/UU2r3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U1.');
fclose(fid);
fid=fopen('./nufftQY/VV2r3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V1.');
fclose(fid);
fid=fopen('./nufftQY/UU2i3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U1.');
fclose(fid);
fid=fopen('./nufftQY/VV2i3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V1.');
fclose(fid);
[fk,ier]=nufft2d3(nj,xj,yj,cj,iflag,eps,nk,sk,tk);


fk=fk(:);
%error1 = norm(fhat-fhat1,2)/norm(fhat,2)
error1 = norm(fhat-fk)/norm(fk)
error2 = norm(fhat1-fk)/norm(fk)
