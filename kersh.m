      function prec=kersh(nequ,nterm,ia,ja,sysmat)
% traduzione in MATLAB della subroutine scritta in FORTRAN kersh.f
%
%  OSSERVAZIONI UTILI per il passaggio da un codice FORTRAN ad uno MATLAB
% Nel passare dalle istruzioni FORTRAN a quelle MATLAB occorre tener presente:
% + in FORTRAN l'istruzione  ELSE IF si puo' scrivere come ELSE IF o ELSEIF
%   in MATLAB  l'analoga istruzione deve essere scritta tutta attaccata
%    ELSEIF altrimenti si apre un altro ciclo IF
% + in FORTRAN un ciclo del tipo DO k=i,j ... END DO dichiara la variabile
%   k come variabile intera e, appena entra nel ciclo la pone uguale a i
%   k=i. Se j<i il ciclo non da' risultati in uscita ma la variabile k
%   rimane comunque del valore uguale a i.
%   Invece, quando i<j, il ciclo DO va avanti e k viene incrementato  di
%   una unita'. Quando k=j, si "lavora" per l'ultima volta dentro il ciclo
%   perche' poi k viene ancora incrementato di un'unita' e quindi il 
%   suo valore diventa k=j+1 >j e  si esce dal ciclo. Il valore in
%   uscita di k e' dunque j+1
%   In MATLAB, invece, una volta che si esce da un ciclo do k=i:j
%   la variabile k e' come cancellata e quindi non puo' essere utilizzata
%   come variabile all'esterno del ciclo.
%   Da qui si puo' osservare come la variabile k2, nel programma FORTRAN
%   e' usata all'esterno del ciclo do k2=i,j-1 con il valore k2=j
%   (sia se i=j sia per i<j)
%         if(j.ge.i) then
%            prec(ia(ja(j))) = prec(ia(ja(j))) + prec(k2)**2
%        end if
%   Lo stesso discorso non si puo' fare in MATLAB perche' k2 viene cancellato
%   al di fuori del ciclo ma le istruzioni analoghe diventano
%         if  j >= i
%              prec(ia(ja(j))) = prec(ia(ja(j))) + prec(j)^2
%        end
%
%
%
%
%  nequ   : numero di equazioni del sistema (n)
%  nterm  : numero di coefficienti non nulli della matrice (nt)
%  ia     : vettore topologico 
%  ja     : vettore con gli indici di colonna
%  sysmat : vettore per la memorizzazione compatta della matrice del sistema
%
%
      zero=0.0;

      for k=1:nterm
         prec(k) = zero;
      end 

      for kk=1:nequ-1;
         k = ia(kk);
         a = sysmat(k) - prec(k);
         if a<=zero
            display('attenzione')
             kk,k,a
             prec(ia(kk-1))
            a = (prec(ia(kk-1)))^2
         end 
         prec(k) = sqrt(a);
 
         i = ia(kk) + 1;
         j = ia(kk+1) - 1;


         for k1 = i:j;
            prec(k1) = (sysmat(k1)-prec(k1))/prec(k);
         end 

         for k2 = i:j-1;
            j1 = ia(ja(k2));
            prec(j1) = prec(j1) + prec(k2)^2;
            i1 = k2 + 1;
            j1 = j1 + 1;
            while (j1<ia(ja(k2)+1) & i1<=j);
               if(ja(j1)==ja(i1));
                  prec(j1) = prec(j1) + prec(k2)*prec(i1);
                  i1 = i1 + 1;
                  j1 = j1 + 1;
               elseif (ja(j1)<ja(i1));
                  j1 = j1 + 1;
               elseif (ja(j1)>ja(i1));
                   i1 = i1 + 1;
               end
            end 
         end 

         if  j >= i  ;
               prec(ia(ja(j))) = prec(ia(ja(j))) + prec(j)^2;
         end 

      end 

      k = ia(nequ);
      a = sysmat(k)-prec(k);
      if(a<=zero) 
        display('attenzione')
          nequ ,a
          prec(ia(nequ-1))
         a = (prec(ia(nequ-1)))^2
      end 
      prec(k) = sqrt(a);

