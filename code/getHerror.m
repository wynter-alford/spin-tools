function Herror = getHerror(testUnitary,targetUnitary,dim,tCyc)
    Uerror = testUnitary*targetUnitary';
    [vecs,valMat] = eig(Uerror);
    
    Herror = zeros(dim,dim);
    for k=1:dim
        ek = (-1i/tCyc)*log(valMat(k,k));
        Herror = Herror + ek*(vecs(:,k)*vecs(:,k)');
    end
end