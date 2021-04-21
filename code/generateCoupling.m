function as = generateCoupling(Nin)
a1 = abs(randn(Nin));
a1 = triu(a1,1) + triu(a1,1)';
a2 = abs(randn(Nin));
a2 = triu(a2,1) + triu(a2,1)';
a3 = abs(randn(Nin));
a3 = triu(a3,1) + triu(a3,1)';
a4 = abs(randn(Nin));
a4 = triu(a4,1) + triu(a4,1)';
a5 = abs(randn(Nin));
a5 = triu(a5,1) + triu(a5,1)';
a6 = abs(randn(Nin));
a6 = triu(a6,1) + triu(a6,1)';
a7 = abs(randn(Nin));
a7 = triu(a7,1) + triu(a7,1)';
a8 = abs(randn(Nin));
a8 = triu(a8,1) + triu(a8,1)';

as={a1,a2,a3,a4,a5,a6,a7, a8};

end