addpath(genpath(pwd))
format long
a = 128;
n1=a;
n2=a;
ms=a;
mt=a;
nj=n1*n2;
eps=1e-12;
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
%xj=pi*rand(nj,1);
%yj=pi*rand(nj,1);
%cj=pi*rand(nj,1);
x=[xj(:) yj(:)];
[k1,k2] = ndgrid((-ms/2):(ms/2-1));
k=[k1(:) k2(:)];


fftconst = iflag*1i/ms*2*pi;
ratiofun = @(k,x)exp(fftconst*k*(x-round(x))');
[U,V] = lowrank(k,x/2/pi*ms,ratiofun,eps,100,100);
xsub = mod(round(x/2/pi*ms),ms)+1;
xxsub = sub2ind([ms ms],xsub(:,1),xsub(:,2));
spPerm = sparse(xxsub,1:nj,ones(1,nj),ms^2,nj);
r = size(V,2)
[n,ncol] = size(cj);
M = repmat(conj(V),1,ncol).*reshape(repmat(cj,r,1),n,r*ncol);
MM = reshape(spPerm*M,ms,ms,r*ncol);
MMM = fft2(MM);
MMMM=fftshift(fftshift(MMM,1),2);
fhat = sum(U.*reshape(MMMM,ms*ms,r),2);
fhat = fhat/nj;

[fk,ier]=nufft2d1(nj,x(:,1),x(:,2),cj,iflag,eps,ms,ms);


fk=fk(:);
%error1 = norm(fhat-fhat1,2)/norm(fhat,2)
error2 = norm(fhat-fk)/norm(fk)

fid=fopen('./nufftQY/U2r1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U.');
fclose(fid);
fid=fopen('./nufftQY/V2r1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V.');
fclose(fid);
fid=fopen('./nufftQY/U2i1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U.');
fclose(fid);
fid=fopen('./nufftQY/V2i1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V.');
fclose(fid);

