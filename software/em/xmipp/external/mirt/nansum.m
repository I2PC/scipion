function y=nansum(x)

x(isnan(x))=0;
y=sum(x);
