FloquetTerms = {};
FloquetUnitaries = {};
for floqTermInd = 0:maxFloqN
    FloquetTerms{floqTermInd+1}=HFloquet(floqTermInd,Taus,toggledHsys,dim); %#ok<*SAGROW> 
    FloquetUnitaries{floqTermInd+1}=expm(1i*FloquetTerms{floqTermInd+1}*2*pi*tCyc);
    raw_QOCT(paramValInd,coupMatInd,floqTermInd+1) = overlap(FloquetUnitaries{floqTermInd+1},FloquetUnitaries{1},N);
    raw_QDCT(paramValInd,coupMatInd,floqTermInd+1) = uDist(FloquetUnitaries{floqTermInd+1},FloquetUnitaries{1},N);

    raw_QOFT(paramValInd,coupMatInd,floqTermInd+1) = overlap(expUnitary,FloquetUnitaries{floqTermInd+1},N);
    raw_QOIT(paramValInd,coupMatInd,floqTermInd+1) = overlap(deltaUnitary,FloquetUnitaries{floqTermInd+1},N);
    raw_QDFT(paramValInd,coupMatInd,floqTermInd+1) = uDist(expUnitary,FloquetUnitaries{floqTermInd+1},N);
    raw_QDIT(paramValInd,coupMatInd,floqTermInd+1) = uDist(deltaUnitary,FloquetUnitaries{floqTermInd+1},N);

    if pulseDivs > 1
        floquetTermsP{floqTermInd+1}=UFloquet(floqTermInd,Taus,toggledHsysP,dim);
        raw_QOCP(paramValInd,coupMatInd,floqTermInd+1) = overlap(floquetTermsP{floqTermInd+1},FloquetUnitaries{1},N);
        raw_QDCP(paramValInd,coupMatInd,floqTermInd+1) = uDist(floquetTermsP{floqTermInd+1},FloquetUnitaries{1},N);

        raw_QOFP(paramValInd,coupMatInd,floqTermInd+1) = overlap(expUnitary,floquetTermsP{floqTermInd+1},N);
        raw_QOIP(paramValInd,coupMatInd,floqTermInd+1) = overlap(deltaUnitary,floquetTermsP{floqTermInd+1},N);
        raw_QDFP(paramValInd,coupMatInd,floqTermInd+1) = uDist(expUnitary,floquetTermsP{floqTermInd+1},N);
        raw_QDIP(paramValInd,coupMatInd,floqTermInd+1) = uDist(deltaUnitary,floquetTermsP{floqTermInd+1},N);
    end
end