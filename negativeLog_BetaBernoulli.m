function f = negativeLog_BetaBernoulli(k,x,a)
    f = 0;
    i = []; 
    if length(a) == length(x)
        for i = 1:length(a)
            k = k;
            %k = k/(k-1);
            f1 = -log(betapdf(x(i),k*a(i)+1,k*(25+1-a(i))));
            f = f + f1;
        end
    else
        error('n° sequences do not match')
    end
end