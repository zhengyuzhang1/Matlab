function [C, Ceq] = nonlcon(x, maxtotaltime)

Ceq = [];
C = (2 * sum(x(2:2:end-1)))*x(end) - maxtotaltime;
end