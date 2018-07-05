function [rv]=nbody(dt,rvmg)
%r0_and_v0=[r0,v0]
l=(length(rvmg)-1)/7;
r0=rvmg(1:l*3); %r0=[m_i,x,y...z]
r10=rvmg(l*3+1:l*6); %v0=[m_i,dx/dt,dy/dt...dz/dt]
m=rvmg(l*6+1:l*7);
G=rvmg(end);
%r2=zeros(length(m),size(r0,2));
for i=1:l %object number
    ii=3*[1:l]-3+i;
    m_ni=m(1:end~=i); %masses other m_i
    r_ni=r0; r_ni(ii)=[]; %radi other than r_i
    r_i=r0(ii); %current postion mass subject to a force
    d=r_i-r_ni; %distnaces from other masses
    a=-G*m_ni./sum(abs(d').^3).*d'; %force/m_i
    r2(i,:)=sum(a,2); %second derivative of r in the j dirrection
end
r1=r10+r2*dt; %r1=r1+dr1=r1+r2*dt
r=r0+r10*dt+.5*r2*dt^2; %r_new=r+dr=r+r1_new*dt+.5*r2*dt^2 x=x+vt+.5at^2
rv=[r,r1,m,G];
end