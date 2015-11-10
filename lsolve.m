      function pvec=lsolve(nequ,nterm,ia,ja,lmat,vec)
%
%   pvec:=(LL^T)^-1*vec (solves the system LL^T*pvec=vec)
%

      zero=0.0;

      for k = 1:nequ;
         pvec(k) = zero;
      end 
      for k = 1:nequ;
         i = ia(k);
         j = ia(k+1) - 1;
         pvec(k) = (vec(k) - pvec(k))/lmat(i);
         for m = i+1:j;
            pvec(ja(m)) = pvec(ja(m)) + lmat(m)*pvec(k);
         end 
      end 
      for k = 1:nequ ;
         n1 = nequ-k+1;
         a = zero;
         i = ia(n1);
         j = ia(n1+1) - 1;
         for m = i:j-1;
            mm = j-m+i;
            a = a + lmat(mm)*pvec(ja(mm));
         end 
         pvec(n1) = (pvec(n1) - a)/lmat(i);
      end

