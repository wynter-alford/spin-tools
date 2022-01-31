% which phase transient to associate with which pulse 
% (note that alpha canbe negative)
function trP = pulseError(pulseIn,roterror,trans,X,Y)  

    if pulseIn == X
        trP = roterror*pulseIn+trans*Y;
    elseif pulseIn == Y
        trP = roterror*pulseIn+trans*(-X);
    elseif pulseIn == -X
        trP = roterror*pulseIn+trans*(-Y);
    elseif pulseIn == -Y
        trP = roterror*pulseIn+trans*X;
    else
        trP = roterror*pulseIn; % transients not yet implemented for pi pulses
    end
end