function y=nanmean(x)

p=isnan(x);
x(p)=0;
y=sum(x(:))/sum(~p(:));
