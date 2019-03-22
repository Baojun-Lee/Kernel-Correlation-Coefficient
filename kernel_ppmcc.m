function R=kernel_ppmcc(x,y)
   Len=length(x);
   x1=sort(x);
   y1=sort(y);
   ir1=x1(floor(Len*0.75))-x1(floor(Len*0.25));
   ir2=y1(floor(Len*0.75))-y1(floor(Len*0.25));

   kernel_size1=5*ir1; 
   kernel_size2=5*ir2; 
  
  % kernel_size1=25*min(std(x),ir1)*Len^(-0.2); %Silverman's rule
  % kernel_size2=25*min(std(y),ir2)*Len^(-0.2); %Silverman's rule
   
   x_median=mean(x1(floor(Len*0.4):floor(Len*0.6)));
   y_median=mean(y1(floor(Len*0.4):floor(Len*0.6)));

   x_ker=gauss_kernel(x,x_median,kernel_size1);
   y_ker=gauss_kernel(y,y_median,kernel_size2);

   x1=(x-x_median).*x_ker;
   y1=(y-y_median).*y_ker;
    
   R=x1'*y1/sqrt(sum(x1.^2)*sum(y1.^2));
   if isnan(R) || isinf(R)
      R=0;
   end
        
        