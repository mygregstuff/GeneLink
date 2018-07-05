function eqn=Lag(L,x,l)
T=diff(diff(L,diff(x,l)),l)
V=diff(l,x)
eqn=0==T-V
end