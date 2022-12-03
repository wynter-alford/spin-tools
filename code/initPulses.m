if phaseTrans ~= 0
    xtransient=expm(-1i*X*pi/2*phaseTrans);
    ytransient=expm(-1i*Y*pi/2*phaseTrans);
    xbtransient=expm(1i*X*pi/2*phaseTrans);
    ybtransient=expm(1i*Y*pi/2*phaseTrans);
else
    xtransient=1;
    ytransient=1;
    xbtransient=1;
    ybtransient=1;
end

UDx=ytransient*expm(-1i*X*rotationError*pi/2)*ytransient;
UDy=xbtransient*expm(-1i*Y*rotationError*pi/2)*xbtransient;
UDxbar=ybtransient*expm(1i*X*rotationError*pi/2)*ybtransient;
UDybar=xtransient*expm(1i*Y*rotationError*pi/2)*xtransient;

Ux=ytransient*expm(-1i*2*pi*(Hsys+f1*X*rotationError)*pulse)*ytransient;
Uy=xbtransient*expm(-1i*2*pi*(Hsys+f1*Y*rotationError)*pulse)*xbtransient;
Uxbar=ybtransient*expm(1i*2*pi*(Hsys+f1*X*rotationError)*pulse)*ybtransient;
Uybar=xtransient*expm(1i*2*pi*(Hsys+f1*Y*rotationError)*pulse)*xtransient;

UDtau = expm(-1i*Hsys*tau*2*pi);
UD2tau = UDtau*UDtau;

Utau = expm(-1i*Hsys*(tau-pulse)*2*pi);
U2tau = expm(-1i*Hsys*(2*tau-pulse)*2*pi);

UdivX=expm(-1i*X*pi/(2*pulseDivs));
UdivY=expm(-1i*Y*pi/(2*pulseDivs));
UdivXbar=expm(1i*X*pi/(2*pulseDivs));
UdivYbar=expm(1i*Y*pi/(2*pulseDivs));