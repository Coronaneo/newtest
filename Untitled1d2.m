addpath(genpath(pwd))
format long
a=16384;
nj=a;
ms=a;
xj=(1:nj)'*pi/8;
cj=exp(1i*(1:nj)/nj)';
k=0:(ms-1);
k=k';
k1=-ms/2:(ms/2-1);
k1=k1';
tol=1e-12;
fhat1=DeCom_NUFFT1D_II(cj,xj/nj,k1,tol);
fhat2=nufft1d2(nj,xj/nj*2*pi,-1,tol,ms,cj);
error=norm(fhat1-fhat2)/norm(fhat1)
fftconst = -1*1i/ms*2*pi;
fun = @(x,k)exp(fftconst*(x-round(x))*k');
%K=13;

[U,V] = lowrank(xj,k1,fun,tol,15,15);
xsub = mod(round(xj),ms)+1;
r=size(V,2)
Id = sparse(1:nj,xsub,ones(1,nj),nj,nj);
MMM=conj(V).*repmat(cj,[1,r]);
MM=fftshift(MMM,1);
fhat =  sum(U.*(Id*fft(MM, [], 1)),2);
%R=conj(V);
%r=size(V,2);
error=norm(fhat-fhat2)/norm(fhat2)
fid=fopen('./nufftQY/Ur2.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U.');
fclose(fid);
fid=fopen('./nufftQY/Vr2.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V.');
fclose(fid);
fid=fopen('./nufftQY/Ui2.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U.');
fclose(fid);
fid=fopen('./nufftQY/Vi2.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V.');
fclose(fid);
