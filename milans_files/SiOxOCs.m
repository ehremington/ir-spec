n1=1; n12=n1^2; d=1e-9;  L=1; No=4;
angle=75*pi/180; s=sin(angle); c=cos(angle);
polar=1;
s2=s*s;
Z=SiOxRI; Q=exspa; Coeff=coeff;
Mm=size(Z); N=Mm(1);
clear k X KK
hold on
for p=1:N-1
    K=Z(p+1,1); k(p)=K;
end

% for pp=1:11
%     x=1+(pp-1)*0.1;
%     for p=1:N-1
%         K=k(p); 
%         Pol=[Coeff(p,1) Coeff(p,2) Coeff(p,3)];
%         Nf(p)=polyval(Pol,x);
%     end
%     A=real(Nf);
%     plot(k,A);
% end

for p=1:N-1
    c(p)=1;
end
plot(k,c,'k');
