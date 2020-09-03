function [y,xnew]=mohinh(xold,u)
xnew=[1.8 -0.7700;1 0]*xold+[1;0]*u;
y=[1 -0.5]*xold;
end