% ASEN 3128 Lab 3
% Matthew Pabin


function motor_forces = ComputeMotorForces(Zc,Lc,Mc,Nc,R,km)

c = inv([-1,-1,-1,-1;-R/sqrt(2),-R/sqrt(2),R/sqrt(2),R/sqrt(2); ...
    R/sqrt(2),-R/sqrt(2),-R/sqrt(2),R/sqrt(2); km,-km,km,-km]);

motor_forces = c .* [Zc; Lc; Mc; Nc];

end

