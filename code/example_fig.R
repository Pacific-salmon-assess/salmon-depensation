Sx=seq(0,1000)
c=3
p_rec3 = ((k*(Sx/Sk)^c)*exp(c*(1-Sx/Sk)))
c=5
p_rec5 = ((k*(Sx/Sk)^c)*exp(c*(1-Sx/Sk)))
c=1.5
p_rec1.5 = ((k*(Sx/Sk)^c)*exp(c*(1-Sx/Sk)))
c=1
p_rec1 = ((k*(Sx/Sk)^c)*exp(c*(1-Sx/Sk)))
c=0.5
p_rec0.5 = ((k*(Sx/Sk)^c)*exp(c*(1-Sx/Sk)))

plot(S, exp_rec,type='n',bty='l')
lines(p_rec0.5~Sx)
lines(p_rec1~Sx)
lines(p_rec1.5~Sx)
lines(p_rec3~Sx)
lines(p_rec5~Sx)
text(100,((k*(100/Sk)^1)*exp(1*(1-100/Sk))),'c=1')
text(100,((k*(100/Sk)^0.5)*exp(0.5*(1-100/Sk))),'c=0.5')
text(100,((k*(100/Sk)^1.5)*exp(1.5*(1-100/Sk))),'c=1.5')
text(100,((k*(100/Sk)^3)*exp(3*(1-100/Sk))),'c=3')
text(100,((k*(100/Sk)^5)*exp(5*(1-100/Sk))),'c=5')