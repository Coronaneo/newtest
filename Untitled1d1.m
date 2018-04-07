addpath(genpath(pwd))
format long
a=131072;
nj=a;
ms=a;
xj=(1:nj)'*pi/8;
cj=exp(1i*(1:nj)/nj)';
k=0:(ms-1);
k=k';
k1=-ms/2:(ms/2-1);
k1=k1';
tol=1e-12;
fhat1=DeCom_NUFFT1D_I(cj,xj/nj,k1,tol);
fhat2=nufft1d1(nj,xj/nj*2*pi,cj,-1,tol,ms);
error=norm(fhat1-fhat2*nj)/norm(fhat1)
fftconst = -1*1i/ms*2*pi;
fun = @(k,x)exp(fftconst*k*(x-round(x))');
%K=13;

[U,V] = lowrank(k1,xj,fun,tol,50,50);
xsub = mod(round(xj),ms)+1;
r=size(V,2)
Id = sparse(xsub,1:nj,ones(1,nj),ms,nj);
MMM=conj(V).*repmat(cj,[1,r]);
MM=Id*(MMM);
M=fftshift(fft(MM, [], 1),1);

fhat =  sum(U.*M,2);
%R=conj(V);
%r=size(V,2);
error1=norm(fhat-fhat2*nj)/norm(fhat)

fid=fopen('./nufftQY/Ur1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U.');
fclose(fid);
fid=fopen('./nufftQY/Vr1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V.');
fclose(fid);
fid=fopen('./nufftQY/Ui1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U.');
fclose(fid);
fid=fopen('./nufftQY/Vi1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V.');
fclose(fid);
