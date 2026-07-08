x=1.15; td=10; No=2;
n1=1; n12=n1^2; d=1e-9;  L=1; 
angle=80*pi/180; s=sin(angle); c=cos(angle);
polar=1; td=d*td;
s2=s*s;
Z=SiOxRI; Q=exspa; Coeff=coeff;
Mm=size(Z); N=Mm(1);
clear k
hold off
for p=1:N-1
    K=Z(p+1,1); k(p)=K;
    Kau(p)=0.95*14000/K^.75-10;
    Nau(p)=20*Kau(p)/(24+0.1*K);
    NAU(p)=Nau(p)+1i*Kau(p);
end

Mexp=size(Q); Nexp=Mexp(1); temp=Q(1660,No);
for p=1:Nexp
    q(p)=Q(p,1);
    Spa(p)=Q(p,No)-temp;
end

Ni=n1; Ni2=Ni^2; 
for p=1:N-1
    K=k(p);
    Pol=[Coeff(p,1) Coeff(p,2) Coeff(p,3)];
    Nf=polyval(Pol,x); Nf2=Nf^2; 
    ui=sqrt(Ni2-n12*s2); uf=sqrt(Nf2-n12*s2); 
    if polar==1
        r=-(ui*Nf2-uf*Ni2)/(ui*Nf2+uf*Ni2); 
        t=2*Ni*Nf*uf/(Nf2*ui+Ni2*uf);
    else if polar==0
         r=(ui-uf)/(ui+uf);
         t=2*uf/(ui+uf);
         end
    end
    ti=1/t;
    M=[ti -r*ti
      -r*ti ti];
    for ss=1:L
        Ni=Nf; Ni2=Nf2;
        Nf=polyval(Pol,x);
        Nf2=Nf^2;
        ui=uf; uf=sqrt(Nf2-n12*s2); 
        D=exp(2*pi*1i*K*td*ui);
        Di=1/D;
        if polar==1
            r=-(ui*Nf2-uf*Ni2)/(ui*Nf2+uf*Ni2); 
            t=2*Ni*Nf*uf/(Nf2*ui+Ni2*uf);
        else if polar==0
            r=(ui-uf)/(ui+uf);
            t=2*uf/(ui+uf);
            end
        end
        ti=1/t;
        Q=[D*ti -r*ti*Di
         -r*ti*D Di*ti];
        M=Q*M;
    end
    Ni=Nf; Ni2=Nf2;
    Nf=NAU(p); 
    Nf2=Nf^2; 
     ui=uf; uf=sqrt(Nf2-n12*s2); 
    D=exp(2*pi*1i*K*td*ui); Di=1/D;
    if polar==1
        r=-(ui*Nf2-uf*Ni2)/(ui*Nf2+uf*Ni2); 
        t=2*Ni*Nf*uf/(Nf2*ui+Ni2*uf);
    else if polar==0
        r=(ui-uf)/(ui+uf);
        t=2*uf/(ui+uf);
        end
    end
    ti=1/t;
    Q=[D*ti -r*ti*Di
     -r*ti*D Di*ti];
    M=Q*M;
    Ni=n1; Ni2=n12; ui=sqrt(n12-n12*s2);
    if polar==1
        r=-(ui*Nf2-uf*Ni2)/(ui*Nf2+uf*Ni2); 
    else if polar==0
        r=(ui-uf)/(ui+uf);
        end
    end
    r14=r;
    A(p)=-2*log10(abs(M(2,1)/(M(2,2)*r14)));
end
A=A';

    
