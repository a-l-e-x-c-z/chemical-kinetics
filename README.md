# chemical-kinetics
matlab/octave script for fitting experimental data to chemical kinetics equations

## example, set of chemical reactions

```
dAdt = -kAB*A -kAZ*A+ kBA*B;
dBdt = -kBA*B - kBC*B - kBE*B + kAB*A;
dCdt = -kCD*C + kBC*B;
dDdt = kCD*C;
dEdt = -kEF*E + kBE*B;
dFdt = kEF*E;
dZdt = kAZ*A;
```
