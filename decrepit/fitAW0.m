function fp = fitAW0(testVars,results)
    %aW0 = fittype('b-(W/W0)^a','independent','W');
    aW0 = fittype('min(1,1+k*(W0-W))','independent','W');
    %aW0 = fittype('min(x+b,3)','independent','x');
    fp = fit(testVars,results,aW0,'StartPoint',[max(testVars)/2,max(testVars)*0.01])
end
