       SUBROUTINE derivs(y,dy,n,a,b)
! Y=Y_i, DY= Y'_i, n-dim
      INTEGER n,k1,k2,k3
      DOUBLE PRECISION y(n),dy(n),a(n,n,n),b(n,n)
      DOUBLE PRECISION S,R,BE,AA
      AA=1.D0
      dy=0.d0
      do k1=1,n
         dy(k1)=-y(k1)**3
!         dy(k1)=-y(k1)*dabs(y(k1))*n*n
          do k2=1,n
            dy(k1)=dy(k1) + b(k1,k2)*y(k2) 
            do k3=1,n
            dy(k1)=dy(k1) + a(k1,k2,k3)*y(k2)*y(k3)  
            end do   
            end do
            end do


      RETURN
       END SUBROUTINE DERIVS
         
