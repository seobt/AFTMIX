function fun = test_linear (b,x)
p=length(b);
%fun=b(1)+b(2)*x(:,1)+b(3)*x(:,2)+b(4)*x(:,3)+b(5)*x(:,4)+b(6)*x(:,5)+b(7)*x(:,6);
fun=b(1)+x(:,[1:p-1])*b(2:p);
