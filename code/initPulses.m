if phaseTrans ~= 0
    xtrans=expm(-1i*X*pi/2*phaseTrans);
    ytrans=expm(-1i*Y*pi/2*phaseTrans);
    xbtrans=expm(1i*X*pi/2*phaseTrans);
    ybtrans=expm(1i*Y*pi/2*phaseTrans);
else
    xtrans=1;
    ytrans=1;
    xbtrans=1;
    ybtrans=1;
end

UDx=ytrans*expm(-1i*X*rotationError*pi/2)*ytrans;
UDy=xbtrans*expm(-1i*Y*rotationError*pi/2)*xbtrans;
UDxbar=ybtrans*expm(1i*X*rotationError*pi/2)*ybtrans;
UDybar=xtrans*expm(1i*Y*rotationError*pi/2)*xtrans;

Ux=ytrans*expm(-1i*2*pi*(Hsys+f1*X*rotationError)*pulse)*ytrans;
Uy=xbtrans*expm(-1i*2*pi*(Hsys+f1*Y*rotationError)*pulse)*xbtrans;
Uxbar=ybtrans*expm(1i*2*pi*(Hsys+f1*X*rotationError)*pulse)*ybtrans;
Uybar=xtrans*expm(1i*2*pi*(Hsys+f1*Y*rotationError)*pulse)*xtrans;

UDtau = expm(-1i*Hsys*tau*2*pi);
UD2tau = UDtau*UDtau;

Utau = expm(-1i*Hsys*(tau-pulse)*2*pi);
U2tau = expm(-1i*Hsys*(2*tau-pulse)*2*pi);
