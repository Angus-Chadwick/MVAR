function [out,dif]=isequalwithequalnans_tolerance(data1,data2,tol)

bothnan=isnan(data1)&isnan(data2);

data1(bothnan)=0;
data2(bothnan)=0;

dif=abs(data1-data2)>tol;

out=not( any(dif(:)') );

