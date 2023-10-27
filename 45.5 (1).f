          double precision function beta(T)
          double precision , intent(in) :: T
          double precision :: tc
          double precision , dimension(100) :: b
          integer:: i
          b(1)=-0.6E-4
          b(2)=0.1E-4
          b(3)=0.9E-4
          b(4)=1.5E-4
          b(5)=2.1E-4
          b(6)=2.6E-4
          b(7)=3.0E-4
          b(8)=3.4E-4
          b(9)=3.8E-4
          b(10)=4.15E-4
          b(11)=4.5E-4
          b(12)=4.8E-4
          b(13)=5.1E-4
          b(14)=5.4E-4
          b(15)=5.7E-4
          b(16)=5.95E-4
          b(17)=6.2E-4
          b(18)=6.45E-4
          b(19)=6.7E-4
          b(20)=6.9E-4
          b(21)=7.1E-4
          tc=T-273.d0
          do i=1,21
            if((tc>=5.d0*(i-1)).and.(tc<5.d0*i))then
         beta=b(i)+((tc-5.d0*real(i-1))/5.d0)*(b(i+1)-b(i))
            endif

          end do
          
          end function beta
          
          double precision function meubl(T)
          double precision , intent(in) :: T
          double precision :: C0,C1,C2,C3,C4
          C0=0.437721d0
          C1=-5.0587E-3
          C2=2.20399E-5
          C3=-4.28147E-8
          C4=3.12531E-11
          meubl=C0+C1*T+C2*T*T+C3*T**3+C4*T**4
          end function meubl
          
           double precision function robl(T)
          double precision , intent(in) :: T
          double precision :: C0,C1,C2
          C0=759.3d0
          C1=1.8475
          C2=-3.52615E-3
          robl=C0+C1*T+C2*T*T
          end function robl
          
          double precision function knf(phi,T)
          double precision , intent(in) :: phi,T
          double precision :: kbf,C0,C1,C2
          C0=-0.748304d0
          C1=7.42386E-3
          C2=-9.64742E-6
          kbf=C0+C1*T+C2*T*T

          knf=kbf*(1+7.47d0*phi)
          end function knf

          double precision function ronf(phi,T)
          double precision , intent(in) :: phi,T
          double precision :: ronp, robf,C0,C1,C2
          ronp=3880.d0
          C0=759.3d0
          C1=1.8475d0
          C2=-3.52615E-3
          robf=C0+C1*T+C2*T*T

          ronf=ronp*phi+(1.d0-phi)*robf
          end function ronf

          double precision function meunf(phi,T)
          double precision , intent(in) :: phi,T
          double precision :: meubf,C0,C1,C2,C3,C4
          C0=0.437721d0
          C1=-5.0587E-3
          C2=2.20399E-5
          C3=-4.28147E-8
          C4=3.12531E-11
          meubf=C0+C1*T+C2*T*T+C3*T**3+C4*T**4
          meunf=meubf*(1+39.11d0*phi+533.9d0*phi*phi)
          end function meunf

          double precision function cnf(phi,T)
          double precision , intent(in) :: phi,T
          double precision :: ronp,robf,cnp,cbf,C0,C1,C2
          ronp=3880.d0
          robf=759.3d0+1.8475d0*T-(3.52615E-3)*T*T
          cnp=773.d0
          C0=5.52022E3
          C1=-8.42432d0
          C2=1.32352E-2
          cbf=C0+C1*T+C2*T*T
          cnf=(phi*ronp*cnp+(1-phi)*robf*cbf)/(phi*ronp+(1-phi)*robf)
          end function cnf

          double precision function Dbnf(T)
          double precision , intent(in) ::T
          double precision :: C0,C1,C2,C3,C4,kb,dp,pi,meubf
          pi=3.141592654d0
          kb=1.38065E-23
          dp=13E-9
          C0=0.437721d0
          C1=-5.0587E-3
          C2=2.20399E-5
          C3=-4.28147E-8
          C4=3.12531E-11
          meubf=C0+ C1*T + C2*T*T + C3*T**3 + C4*T**4
          Dbnf=kb*T/(3.d0*pi*meubf*dp)
          end function Dbnf

          double precision function Dtlnf(phi,T)
          double precision , intent(in) ::phi,T
          double precision :: rob,meubf,kbf,C0,C1,C2,C3,C4,
     &    Ck0,Ck1,Ck2,kp
          C0=0.437721d0
          C1=-5.0587E-3
          C2=2.20399E-5
          C3=-4.28147E-8
          C4=3.12531E-11
          meubf=C0 + C1*T + C2*T*T + C3*T**3 + C4*T**4
          robf=759.3d0+1.8475d0*T-(3.52615E-3)*T*T
          Ck0=-0.748304d0
          Ck1=7.42386E-3
          Ck2=-9.64742E-6
          kbf=Ck0 + Ck1*T + Ck2*T*T
          kp=36.d0
          Dtlnf=0.26d0*(kbf/(2.d0*kbf+kp))*(meubf/robf)*(phi/T)
          end function Dtlnf



         program smple_cmp
         implicit none
          double precision ,dimension(0:500,0:500):: u,v,p,pe,theta,
     &    phi,k,ro,
     &   meu,c,phiE,phiW,phiN,phiS,phiP,phiB,uW,uE,uN,uS,uP,uB,vW,vu,uv,
     &   vE,vS,vN,vP,vB,thW,thE,thN,thS,thP,thB,uPP,vPP,rov,rou,
     &   meuu,meuv,phieu,phir,eu,ur,ev,vr,eth,thr,pee,peu,pepx,pepy,pep,
     &   er,epp,phiuu,uuoold,phiold,thold,thetaee,phieee,pepys,pepyn,
     &   vuw,vue,us_grid,vs_grid,Db,Dtl,betabf,meubfl,robfl,
     &   nr,si,vold,robflu,betabfu

         double precision ,dimension(0:500):: dytheta,peul,pepl,peub,
     &    pepb,hl,ys_grid,Nul,ddy,ReL,RiL,K1,GrL

         double precision ::L,d,dx,dy,dt,U0,meu0,Re0,Nbt0,Db0,Dtl0,Sc0,
     &   Le0,Pr0,T0,Ts,phi0,ro0,theta0,dp,bitabf,g,Reb,rop,kp,cp,c0,
     &   robf,meubf,kbf,k0,Grbf,F,Patm,delT,Ri,r,phipe,
     &   phiu,upe,uiu,vpe,viu,thpe,thiu,pi,meunf,knf,cnf,
     &   ronf,peuc,pepc,e,err,ep,A,kb,ueebig,phibig,thbig,
     &   dterm,nondterm,dtt,wup,pup,eup,svp,pvp,nvp,wuv,puv,euv,svv,
     &   pvv,nvv,wu,pu,euup,sv,pv,nvvp,phid,havg,sum,rp,rv1,rv2,dy1,dy2,
     &  Nuavg,tolerance,us_sum,phis_sum,thetas_sum,u_rmse,phi_rmse,
     &  theta_rmse,dyp,dy3,robitanf,bitanp,Rinf,Grnf,Dbnf,Dtlnf,Tl,
     &  beta,robl,meubl,veebig,vs_sum,v_rmse,tolerance1,Ribf

         integer:: N,M,i,j,ki,t,flag,checkn,us_num,
     &  thetas_num,phis_num,M1,M2,M3,MZ,vs_num
         !! properties values
         pi=3.141592654
         kb=1.38065E-23
         dp=13e-9
         Ts=368.d0
         T0=298.d0
         phi0=0.05d0
         
         rop=3880.d0
         kp=36.d0
         cp=773.d0
         robf=robl(T0)
         meubf=meubl(T0)
         !kbf=0.6d0
         !bitabf=0.00026d0
         bitanp=0.850E-5
         g=9.806d0
         Patm=101325.d0
         L=1.d0
         d=0.05d0
         
         U0=0.05d0
         !U0=0.05d0

         !! properties
         bitabf=beta(T0)
         meu0=meunf(phi0,T0)
         ro0=ronf(phi0,T0)
         k0=knf(phi0,T0)
         c0=cnf(phi0,T0)

         robitanf=(robf*bitabf)*(1.d0-phi0)+(rop*bitanp)*phi0
         Grnf=(robitanf/ro0)*g*(Ts-T0)*(ro0**2.d0)*(L**3.d0)/(meu0**2)
         !! reference non dimensional numbers
         Grbf=bitabf*g*(Ts-T0)*(robf**2.d0)*(L**3.d0)/(meubf**2)

         Reb=U0*L*robf/meubf

         Ribf=(Grbf/(Reb**2))

         Re0=U0*L*ro0/meu0
         Rinf=Grnf/(Re0*Re0)
         
         delT=Ts-T0
         Db0=Dbnf(T0)
         Dtl0=Dtlnf(phi0,T0)
         
         Sc0=meu0/(Db0*ro0)
         Le0=k0/(rop*cp*Db0*phi0)
         
         Ri=(Grbf/((Reb)**2.d0))*(1.d0-phi0)
         
         F=real((g*L)/(U0*U0))
         
         Nbt0=phi0*Db0/(Dtl0*(Ts-T0))
         Pr0=meu0*c0/k0


          N=51
         MZ=10
         M1=5
         M2=15
         M3=20
         M=(MZ-1)*M1+M2+M3



            !!!main iteration starts***************************************
            
           ki=0
           
            checkn=0
            !!! main loop starts ===============================================================
            do
               !!!!!!!!!!!!!!!!!!!! initialization

           do i=0,N
            do j=0,M
            u(i,j)=0.d0
            v(i,j)=0.d0

            theta(i,j)=0.d0
            phi(i,j)=1.d0
            p(i,j)=Patm/(ro0*U0*U0)
            pe(i,j)=0.d0

            enddo
           enddo






              dx=1.d0/dfloat(N-1)
         dy1=(d/L)/dfloat(MZ*M2*M3)
         dy2=(d/L)/dfloat(MZ*M2)
         dy3=(d/L)/dfloat(MZ*M1)
         !print*,'dy2',dy2

         !rp=dy2/dy1
            dt=0.001d0*(0.5d0**checkn)
            !dt=0.001d0
            A=dt

            print*, "time step",dt,"attempt",checkn
             ki=0
             tolerance=5.d0
             do while (tolerance>(0.0005*dt))
             !do ki=1,3


            do i=1,N
             do j=1,M
               phid=phi0*phi(i,j)
               Tl=T0+(Ts-T0)*theta(i,j)
               ro(i,j)=ronf(phid,Tl)
               k(i,j)=knf(phid,Tl)
               meu(i,j)=meunf(phid,Tl)
               c(i,j)=cnf(phid,Tl)
               Dtl(i,j)=Dtlnf(phid,Tl)
               Db(i,j)=Dbnf(Tl)
               betabf(i,j)=beta(Tl)
               robfl(i,j)=robl(Tl)
               meubfl(i,j)=meubl(Tl)

             end do
            enddo

            do i=1,N-1
            do j=1,M

            robflu(i,j)=(robfl(i,j)+robfl(i+1,j))/2.d0
            betabfu(i,j)=(betabf(i,j)+betabf(i+1,j))/2.d0

            enddo
            enddo
            
            !!!!! extrapolated values of image points
            
            do j=1,M

               theta(N+1,j)=theta(N-1,j)
               phi(N+1,j)=phi(N-1,j)
               u(N,j)=u(N-1,j)
               v(N+1,j)=v(N-1,j)

             enddo

               do i=1,N
               v(i,M)=v(i,M-1)
               enddo



            !!! calculating coefficients of equation

            !!!!! near slot side u equation

             i=1
          do j=2,M-1

           if (j<(M3+1))then
            rp=1.d0
            dy=dy1
            else if(j==(M3+1))then
            rp=dy2/dy1
            dy=dy1
            else if(j<M2+M3)then
            rp=1.d0
            dy=dy2
            else if(j==M2+M3)then
            dy=dy2
            rp=dy3/dy2
            else
            rp=1.d0
            dy=dy3
           end if                                        !(M/2)                                                   !M-1


          vuw(i,j)=v(i,j-1)+(1.d0/(1.d0+rp))*(v(i,j)-v(i,j-1))
          vue(i,j)=v(i+1,j-1)+(1.d0/(1.d0+rp))*(v(i+1,j)-v(i+1,j-1))

            vu(i,j)=0.5d0*(vuw(i,j)+vue(i,j))
          rou(i,j)= (ro(i,j)+ro(i+1,j))/2.d0
          meuu(i,j)=(meu(i,j)+meu(i+1,j))/2.d0
          meuu(i,j+1)=(meu(i,j+1)+meu(i+1,j+1))/2.d0
          meuu(i,j-1)=(meu(i,j-1)+meu(i+1,j-1))/2.d0
          phiuu(i,j)=(phi(i,j)+phi(i+1,j))/2.d0

         if(u(1,j)>0.d0)then
          wu=-2.d0*u(1,j)
          pu=2.d0*u(1,j)
          euup=0.d0
          else
          wu=0.d0
          pu=-u(1,j)
          euup=u(1,j)
          endif

             if(vu(i,j)>0.d0)then
             sv=-vu(i,j)
             pv=vu(i,j)
             nvvp=0.d0
             else
             sv=0.d0
             pv=-vu(i,j)/rp
             nvvp=vu(i,j)/rp
             endif




          uW(1,j)= wu/dx+(4.d0/(Re0*3.d0*dx*dx*meu0))*(ro0/rou(1,j))*
     & (meu(2,j)-meu(1,j))-(8.d0/(3.d0*Re0*dx*dx*meu0))*(ro0/rou(1,j))*
     &      meuu(1,j)

             uE(1,j)= -(1.d0/(Re0*3.d0*dx*dx*meu0))*(ro0/rou(1,j))*
     &  (meu(2,j)-meu(1,j))-(4.d0/(3.d0*Re0*dx*dx*meu0))*(ro0/rou(1,j))*
     &      meuu(1,j) +euup/dx

       uS(i,j)= (1.d0/(Re0*meu0))*(ro0/rou(i,j))*(meuu(i,j+1)/(rp*dy*
     & (1.d0+rp))-(1.d0-rp*rp)*meuu(i,j)/
     & (rp*dy*(1.d0+rp))-meuu(i,j-1)*rp*rp/(rp*dy*(1.d0+rp)))*rp/
     & (dy*(1.d0+rp))-(1.d0/(Re0*meu0))*(ro0/rou(i,j))*
     &  meuu(i,j)*2.d0/((1.d0+rp)*dy*dy)+sv/dy

       uN(i,j)= nvvp/dy-(1.d0/(Re0*meu0))*(ro0/rou(i,j))*(meuu(i,j+1)/
     & (rp*dy*(1.d0+rp))-(1.d0-rp*rp)*meuu(i,j)/
     &  (rp*dy*(1.d0+rp))-meuu(i,j-1)*
     &  rp*rp/(rp*dy*(1.d0+rp)))/(rp*(dy*(1.d0+rp)))-(1.d0/(Re0*meu0))*
     & (ro0/rou(i,j))*meuu(i,j)*2.d0/(rp*(1.d0+rp)*dy*dy)


       uPP(i,j)=1.d0/dt+pu/dx+pv/dy+(4.d0*meuu(1,j)/(Re0*meu0*dx**2))*
     & (ro0/rou(i,j))-(1.d0/(Re0*dx*dx*meu0))*(ro0/rou(1,j))*
     &  (meu(2,j)-meu(1,j))+(2.d0*meuu(i,j)/(Re0*meu0*rp*dy*dy))*
     &  (ro0/rou(i,j))+(1.d0/(Re0*meu0))*(ro0/rou(i,j))*(meuu(i,j+1)/
     & (rp*dy*(1.d0+rp))-(1.d0-rp*rp)*meuu(i,j)/
     & (rp*dy*(1.d0+rp))-meuu(i,j-1)*
     &  rp*rp/(rp*dy*(1.d0+rp)))*(1.d0-rp*rp)/(rp*dy*(1+rp))

          uB(1,j)=u(1,j)/dt-(ro0/rou(1,j))*(p(2,j)-p(1,j))/dx+
     & Ri*(betabfu(1,j)/bitabf)*(theta(1,j)+theta(2,j))*0.5d0*
     & (robflu(1,j)/rou(1,j))-
     & (F)*phi0*(phiuu(1,j)-1.d0)*(rop-robflu(1,j))/rou(1,j)



            enddo




           !! near the plate v equation


           do i=2,N-1
              j=1
              dy=dy1
              r=(dx/dy)**2
            rov(i,1)=(ro(i,1)+ro(i,2))/2.d0
          meuv(i,1)=(meu(i,1)+meu(i,2))/2.d0
          meuv(i+1,1)=(meu(i+1,1)+meu(i+1,2))/2.d0
          meuv(i-1,1)=(meu(i-1,1)+meu(i-1,2))/2.d0
          uv(i,1)=(u(i,1)+u(i,2)+u(i-1,2)+u(i-1,1))/4.d0

          if(uv(i,1)>0.d0)then
          wu=-uv(i,1)
          pu=uv(i,1)
          euup=0.d0
          else
          wu=0.d0
          pu=-uv(i,1)
          euup=uv(i,1)
          endif

          if(v(i,1)>0.d0)then
          sv=-2.d0*v(i,1)
          pv=2.d0*v(i,1)
          nvvp=0.d0
          else
          sv=0.d0
          pv=-v(i,1)
          nvvp=v(i,1)
          endif


         vW(i,1)= wu/dx+(1.d0/(Re0*4.d0*dx*dx*meu0))*(ro0/rov(i,1))*
     & (meuv(i+1,1)-meuv(i-1,1))-(1.d0/(Re0*dx*dx*meu0))*(ro0/rov(i,1))*
     &  meuv(i,1)

           vE(i,1)= -(1.d0/(Re0*4.d0*dx*dx*meu0))*(ro0/rov(i,1))*
     & (meuv(i+1,1)-meuv(i-1,1))-(1.d0/(Re0*dx*dx*meu0))*(ro0/rov(i,1))*
     &  meuv(i,1) + euup/dx

           vS(i,1)= (4.d0/(Re0*3.d0*dy*dy*meu0))*(ro0/rov(i,1))*
     & (meu(i,2)-meu(i,1))-(8.d0/(3.d0*Re0*dy*dy*meu0))*(ro0/rov(i,1))*
     &  meuv(i,1) + sv/dy

          vN(i,1)= nvvp/dy-(1.d0/(Re0*3.d0*dy*dy*meu0))*(ro0/rov(i,1))*
     &  (meu(i,2)-meu(i,1))-(4.d0/(Re0*3.d0*dy*dy*meu0))*(ro0/rov(i,1))*
     &   meuv(i,1)

          vPP(i,1)=1.d0/dt + pu/dx + pv/dy+(2.d0*meuv(i,1)/(meu0*
     &   Re0*dx**2))*(1.d0+2.d0*r)*(ro0/rov(i,1))-(1/(Re0*dy*dy*meu0))*
     &   (ro0/rov(i,1))*(meu(i,2)-meu(i,1))

          vB(i,1)=v(i,1)/dt-(ro0/rov(i,1))*(p(i,2)-p(i,1))/dy


            enddo

       !  !! near plate  on exit side
       
              i=N
              j=1
              dy=dy1
              r=(dx/dy)**2
            rov(N,1)=(ro(N,1)+ro(N,2))/2.d0
          meuv(N,1)=(meu(N,1)+meu(N,2))/2.d0
          uv(N,1)=(u(N,1)+u(N,2)+u(N-1,2)+u(N-1,1))/4.d0

          if(uv(N,1)>0.d0)then
          wu=-uv(N,1)
          pu=uv(N,1)
          euup=0.d0
          else
          wu=0.d0
          pu=-uv(N,1)
          euup=uv(N,1)
          endif

          if(v(N,1)>0.d0)then
          sv=-2.d0*v(N,1)
          pv=2.d0*v(N,1)
          nvvp=0.d0
          else
          sv=0.d0
          pv=-v(N,1)
          nvvp=v(N,1)
          endif


         vW(N,1)= wu/dx+(1.d0/(Re0*4.d0*dx*dx*meu0))*(ro0/rov(N,1))*
     & (3.d0*meuv(N,1)-4.d0*meuv(N-1,1)+meuv(N-2,1))-
     &  (1.d0/(Re0*dx*dx*meu0))*(ro0/rov(N,1))*meuv(N,1)

           vE(N,1)= -(1.d0/(Re0*4.d0*dx*dx*meu0))*(ro0/rov(N,1))*
     & (3.d0*meuv(N,1)-4.d0*meuv(N-1,1)+meuv(N-2,1))-
     &  (1.d0/(Re0*dx*dx*meu0))*(ro0/rov(N,1))*
     &  meuv(N,1) + euup/dx

           vS(N,1)= (4.d0/(Re0*3.d0*dy*dy*meu0))*(ro0/rov(N,1))*
     & (meu(N,2)-meu(N,1))-(8.d0/(3.d0*Re0*dy*dy*meu0))*(ro0/rov(N,1))*
     &  meuv(N,1) + sv/dy

          vN(N,1)= nvvp/dy-(1.d0/(Re0*3.d0*dy*dy*meu0))*(ro0/rov(N,1))*
     &  (meu(N,2)-meu(N,1))-(4.d0/(Re0*3.d0*dy*dy*meu0))*(ro0/rov(N,1))*
     &   meuv(N,1)

          vPP(N,1)=1.d0/dt + pu/dx + pv/dy+(2.d0*meuv(N,1)/(meu0*
     &   Re0*dx**2))*(1.d0+2.d0*r)*(ro0/rov(N,1))-(1/(Re0*dy*dy*meu0))*
     &   (ro0/rov(N,1))*(meu(N,2)-meu(N,1))

          vB(N,1)=v(N,1)/dt-(ro0/rov(N,1))*(p(N,2)-p(N,1))/dy




            !!! end of near plate on exit v coeff
       



            flag=0
             !! coefficient of phi equation




           !theta(N+1,1)=1.d0
             do i=2,N-1
              j=1
              dytheta(i)=-theta(i,3)+4*theta(i,2)-3*theta(i,1)
       K1(i)=(1.d0/Nbt0)*dytheta(i)*(Dtl(i,1)/Dtl0)*(Db0/Db(i,1))
           phi(i,0)=phi(i,2)+K1(i)

                phiW(i,1)=(-1.d0/(Re0*Sc0*(dx**2)))*(Db(i,1)/Db0)-
     &   1.d0/dx+(1.d0/Db0)*(Db(i+1,1)-Db(i-1,1))/(4.d0*(dx**2.d0)*Re0*
     &   Sc0)

            phiE(i,1)=(-1.d0/(Re0*Sc0*(dx**2)))*(Db(i,1)/Db0)-
     &   (1.d0/Db0)*(Db(i+1,1)-Db(i-1,1))/(4.d0*(dx**2)*Re0*Sc0)

         phiS(i,1)=(-1.d0/(Re0*Sc0))*(Db(i,1)/Db0)*(1.d0/(dy1**2))+
     & (1.d0/(Re0*Sc0*Db0))*(-Db(i,3)+4.d0*Db(i,2)-3.d0*Db(i,1))/
     & (4.d0*dy1**2)

         phiN(i,1)=(-1.d0/(Re0*Sc0))*(Db(i,1)/Db0)*(1.d0/(dy1**2))-
     & (1.d0/(Re0*Sc0*Db0))*(-Db(i,3)+4.d0*Db(i,2)-3.d0*Db(i,1))/
     &  (4.d0*dy1**2)


       phiP(i,1)= 1.d0/dt+1.d0/dx+(2.d0/(Re0*Sc0*dx**2))*(Db(i,1)/Db0)+
     &  (2.d0/(Re0*Sc0*dy1**2))*(Db(i,1)/Db0)
       
     
         phiB(i,1)=phi(i,1)/dt +(1.d0/(Nbt0*Sc0*Re0))*((Dtl(i,1)/Dtl0)*
     &	 ((theta(i+1,1)-2.d0*theta(i,1)+theta(i-1,1))/(dx**2)+
     &   (-theta(i,4)+4.d0*theta(i,3)-5.d0*theta(i,2)+2.d0*theta(i,1))/
     &   (dy1**2))+(1.d0/Dtl0)*((Dtl(i+1,1)-
     &   Dtl(i-1,1))*(theta(i+1,1)-theta(i-1,1))/
     &   (4.d0*dx*dx)+(-Dtl(i,3)+4.d0*Dtl(i,2)-3.d0*Dtl(i,1))*
     &   (-theta(i,3)+4.d0*theta(i,2)-3.d0*theta(i,1))/(4.d0*dy1*dy1)))




             enddo

             !! for corner point
                  i=N
                  j=1
              dytheta(N)=-theta(N,3)+4*theta(N,2)-3*theta(N,1)
        K1(N)=(1.d0/Nbt0)*dytheta(N)*(Dtl(N,1)/Dtl0)*(Db0/Db(N,1))
           phi(N,0)=phi(N,2)+K1(N)

                phiW(N,1)=(-1.d0/(Re0*Sc0*(dx**2)))*(Db(N,1)/Db0)-
     &   1.d0/dx+(1.d0/Db0)*(3.d0*Db(N,1)-4.d0*Db(N-1,1)+Db(N-2,1))/
     &  (4.d0*(dx**2.d0)*Re0*Sc0)

            phiE(N,1)=(-1.d0/(Re0*Sc0*(dx**2)))*(Db(N,1)/Db0)-
     &   (1.d0/Db0)*(3.d0*Db(N,1)-4.d0*Db(N-1,1)+Db(N-2,1))/
     &  (4.d0*(dx**2)*Re0*Sc0)

         phiS(N,1)=(-1.d0/(Re0*Sc0))*(Db(N,1)/Db0)*(1.d0/(dy1**2))+
     & (1.d0/(Re0*Sc0*Db0))*(-Db(N,3)+4.d0*Db(N,2)-3.d0*Db(N,1))/
     & (4.d0*dy1**2)

         phiN(N,1)=(-1.d0/(Re0*Sc0))*(Db(N,1)/Db0)*(1.d0/(dy1**2))-
     & (1.d0/(Re0*Sc0*Db0))*(-Db(N,3)+4.d0*Db(N,2)-3.d0*Db(N,1))/
     &  (4.d0*dy1**2)


       phiP(N,1)= 1.d0/dt+1.d0/dx+(2.d0/(Re0*Sc0*dx**2))*(Db(N,1)/Db0)+
     &  (2.d0/(Re0*Sc0*dy1**2))*(Db(N,1)/Db0)


         phiB(N,1)=phi(N,1)/dt +(1.d0/(Nbt0*Sc0*Re0))*((Dtl(N,1)/Dtl0)*
     &	 ((2.d0*theta(N,1)-5.d0*theta(N-1,1)+4.d0*theta(N-2,1)-
     &   theta(N-3,1))/(dx**2)+
     &   (-theta(N,4)+4.d0*theta(N,3)-5.d0*theta(N,2)+2.d0*theta(N,1))/
     &  (dy1**2))+(1.d0/Dtl0)*((3.d0*Dtl(N,1)-4.d0*Dtl(N-1,1)+
     &  Dtl(N-2,1))*(3.d0*theta(N,1)-4.d0*theta(N-1,1)+theta(N-2,1))/
     &   (4.d0*dx*dx)+(-Dtl(N,3)+4.d0*Dtl(N,2)-3.d0*Dtl(N,1))*
     &   (-theta(N,3)+4.d0*theta(N,2)-3.d0*theta(N,1))/(4.d0*dy1*dy1)))



             !!! corner point end
             
             
            do i=2,N-1
             do j=2,M-1

              if (j<(M3+1))then
            rp=1.d0
            dy=dy1
            else if(j==(M3+1))then
            rp=dy2/dy1
            dy=dy1
            else if(j<M2+M3)then
            rp=1.d0
            dy=dy2
            else if(j==M2+M3)then
            dy=dy2
            rp=dy3/dy2
            else
            rp=1.d0
            dy=dy3
           end if


              up(i,j)=(u(i,j)+u(i-1,j))/2.d0

             vp(i,j)=v(i,j-1)+(1.d0/(1.d0+rp))*(v(i,j)-v(i,j-1))


             if (up(i,j)>0.d0)then
             wup=-up(i,j)
             pup=up(i,j)
             eup=0.d0
             else
             wup=0.d0
             pup=-up(i,j)
             eup=up(i,j)
             endif

             if(vp(i,j)>0.d0)then
             svp=-vp(i,j)
             pvp=vp(i,j)
             nvp=0.d0
             else
             svp=0.d0
             pvp=-vp(i,j)/rp
             nvp=vp(i,j)/rp
             endif

             phiW(i,j)=(-1.d0/(Re0*Sc0*(dx**2)))*(Db(i,j)/Db0)+
     &   wup/dx+(1.d0/Db0)*(Db(i+1,j)-Db(i-1,j))/(4.d0*(dx**2.d0)*Re0*
     &   Sc0)

            phiE(i,j)=(-1.d0/(Re0*Sc0*(dx**2)))*(Db(i,j)/Db0)-
     &   (1/Db0)*(Db(i+1,j)-Db(i-1,j))/(4.d0*(dx**2)*Re0*Sc0)+
     &    eup/dx

         phiS(i,j)=svp/dy+(-1.d0/(Re0*Sc0))*(Db(i,j)/Db0)*
     & (2.d0/((1.d0+rp)*dy**2))+(1.d0/Db0)*(Db(i,j+1)/(rp*(1.d0+rp)*dy)-
     &  Db(i,j)*(1.d0-rp*rp)/(rp*(1.d0+rp)*dy)-rp*Db(i,j-1)/((1.d0+rp)*
     &  dy))*(rp/((1.d0+rp)*dy))/(Re0*Sc0)

          phiN(i,j)=nvp/dy+(-1.d0/(Re0*Sc0))*(Db(i,j)/Db0)*
     & (2.d0/(rp*(1.d0+rp)*dy**2))-(1.d0/Db0)*
     & (Db(i,j+1)/(rp*(1.d0+rp)*dy)-
     &  Db(i,j)*(1.d0-rp*rp)/(rp*(1.d0+rp)*dy)-rp*Db(i,j-1)/
     & ((1.d0+rp)*dy))*(1.d0/(rp*(1.d0+rp)*dy))/(Re0*Sc0)


       phiP(i,j)= 1.d0/dt+pup/dx+pvp/dy+(2.d0/(Re0*Sc0*dx**2))*
     &  (Db(i,j)/Db0)+(2.d0/(rp*Re0*Sc0*dy**2))*(Db(i,j)/Db0)+
     & (1.d0/Db0)*(Db(i,j+1)/(rp*(1.d0+rp)*dy)-Db(i,j)*(1.d0-rp*rp)/
     & (rp*(1.d0+rp)*dy)-rp*Db(i,j-1)/((1.d0+rp)*dy))*(1.d0/(Re0*Sc0))*
     &  (1.d0-rp*rp)/(rp*dy*(1.d0+rp))

         phiB(i,j)=phi(i,j)/dt +(1.d0/(Nbt0*Sc0*Re0))*((Dtl(i,j)/Dtl0)*
     &	 ((theta(i+1,j)-2.d0*theta(i,j)+theta(i-1,j))/(dx**2)+
     &   (theta(i,j+1)*2.d0/(rp*dy*dy*(1.d0+rp))-(2.d0/(rp*dy*dy))*
     &   theta(i,j)+(2.d0/(dy*dy*(1.d0+rp)))*theta(i,j-1)))+(1.d0/Dtl0)*
     & ((Dtl(i+1,j)-Dtl(i-1,j))*(theta(i+1,j)-theta(i-1,j))/
     & (4.d0*dx*dx)+(Dtl(i,j+1)/(rp*(1.d0+rp)*dy)-Dtl(i,j)*(1.d0-rp*rp)/
     &  (rp*(1.d0+rp)*dy)-rp*Dtl(i,j-1)/((1.d0+rp)*dy))*(theta(i,j+1)/
     &  (rp*(1.d0+rp)*dy)-theta(i,j)*(1.d0-rp*rp)/
     &  (rp*(1.d0+rp)*dy)-rp*theta(i,j-1)/((1.d0+rp)*dy))))



        enddo
        enddo

          !!! image point for phi at exit
            i=N
             do j=2,M-1

              if (j<(M3+1))then
            rp=1.d0
            dy=dy1
            else if(j==(M3+1))then
            rp=dy2/dy1
            dy=dy1
            else if(j<M2+M3)then
            rp=1.d0
            dy=dy2
            else if(j==M2+M3)then
            dy=dy2
            rp=dy3/dy2
            else
            rp=1.d0
            dy=dy3
           end if

                                                         
              up(i,j)=(u(i,j)+u(i-1,j))/2.d0

             vp(i,j)=v(i,j-1)+(1.d0/(1.d0+rp))*(v(i,j)-v(i,j-1))


             if (up(i,j)>0.d0)then
             wup=-up(i,j)
             pup=up(i,j)
             eup=0.d0
             else
             wup=0.d0
             pup=-up(i,j)
             eup=up(i,j)
             endif

             if(vp(i,j)>0.d0)then
             svp=-vp(i,j)
             pvp=vp(i,j)
             nvp=0.d0
             else
             svp=0.d0
             pvp=-vp(i,j)/rp
             nvp=vp(i,j)/rp
             endif

             phiW(i,j)=(-1.d0/(Re0*Sc0*(dx**2)))*(Db(i,j)/Db0)+
     &   wup/dx+(1.d0/Db0)*(3.d0*Db(N,j)-4.d0*Db(N-1,j)+Db(N-2,j))/
     &  (4.d0*(dx**2.d0)*Re0*Sc0)

            phiE(i,j)=(-1.d0/(Re0*Sc0*(dx**2)))*(Db(i,j)/Db0)-
     &   (1/Db0)*(3.d0*Db(N,j)-4.d0*Db(N-1,j)+Db(N-2,j))/
     & (4.d0*(dx**2)*Re0*Sc0)+eup/dx

         phiS(i,j)=svp/dy+(-1.d0/(Re0*Sc0))*(Db(i,j)/Db0)*
     & (2.d0/((1.d0+rp)*dy**2))+(1.d0/Db0)*(Db(i,j+1)/(rp*(1.d0+rp)*dy)-
     &  Db(i,j)*(1.d0-rp*rp)/(rp*(1.d0+rp)*dy)-rp*Db(i,j-1)/((1.d0+rp)*
     &  dy))*(rp/((1.d0+rp)*dy))/(Re0*Sc0)

          phiN(i,j)=nvp/dy+(-1.d0/(Re0*Sc0))*(Db(i,j)/Db0)*
     & (2.d0/(rp*(1.d0+rp)*dy**2))-(1.d0/Db0)*
     & (Db(i,j+1)/(rp*(1.d0+rp)*dy)-
     &  Db(i,j)*(1.d0-rp*rp)/(rp*(1.d0+rp)*dy)-rp*Db(i,j-1)/
     & ((1.d0+rp)*dy))*(1.d0/(rp*(1.d0+rp)*dy))/(Re0*Sc0)


       phiP(i,j)= 1.d0/dt+pup/dx+pvp/dy+(2.d0/(Re0*Sc0*dx**2))*
     &  (Db(i,j)/Db0)+(2.d0/(rp*Re0*Sc0*dy**2))*(Db(i,j)/Db0)+
     & (1.d0/Db0)*(Db(i,j+1)/(rp*(1.d0+rp)*dy)-Db(i,j)*(1.d0-rp*rp)/
     & (rp*(1.d0+rp)*dy)-rp*Db(i,j-1)/((1.d0+rp)*dy))*(1.d0/(Re0*Sc0))*
     &  (1.d0-rp*rp)/(rp*dy*(1.d0+rp))

         phiB(i,j)=phi(i,j)/dt +(1.d0/(Nbt0*Sc0*Re0))*((Dtl(i,j)/Dtl0)*
     &	 ((2.d0*theta(N,j)-5.d0*theta(N-1,j)+4.d0*theta(N-2,j)-
     &  theta(N-3,j))/(dx**2)+
     &   (theta(i,j+1)*2.d0/(rp*dy*dy*(1.d0+rp))-(2.d0/(rp*dy*dy))*
     &   theta(i,j)+(2.d0/(dy*dy*(1.d0+rp)))*theta(i,j-1)))+(1.d0/Dtl0)*
     & ((3.d0*Dtl(N,j)-4.d0*Dtl(N-1,j)+Dtl(N-2,j))*
     &  (3.d0*theta(N,j)-4.d0*theta(N-1,j)+theta(N-2,j))/
     & (4.d0*dx*dx)+(Dtl(i,j+1)/(rp*(1.d0+rp)*dy)-Dtl(i,j)*(1.d0-rp*rp)/
     &  (rp*(1.d0+rp)*dy)-rp*Dtl(i,j-1)/((1.d0+rp)*dy))*(theta(i,j+1)/
     &  (rp*(1.d0+rp)*dy)-theta(i,j)*(1.d0-rp*rp)/
     &  (rp*(1.d0+rp)*dy)-rp*theta(i,j-1)/((1.d0+rp)*dy))))



        enddo




         do i=2,N
         do j=1,M-1
         dterm=abs(phiP(i,j))
         nondterm=abs(phiE(i,j))+abs(phiW(i,j))+abs(phiN(i,j))+
     &   abs(phiS(i,j))
         if (dterm<nondterm)then
         flag=1
         print*,"phi dia dom violated"
         endif
         enddo
         enddo

         if(flag>0)exit


         !





          !! coefficient of u equation
               do i=2,N-1
          do j=2,M-1                                                   !M-1

            if (j<(M3+1))then
            rp=1.d0
            dy=dy1
            else if(j==(M3+1))then
            rp=dy2/dy1
            dy=dy1
            else if(j<M2+M3)then
            rp=1.d0
            dy=dy2
            else if(j==M2+M3)then
            dy=dy2
            rp=dy3/dy2
            else
            rp=1.d0
            dy=dy3
           end if

           !j=(M/2)+1

          !dy=dy1
          vuw(i,j)=v(i,j-1)+(1.d0/(1.d0+rp))*(v(i,j)-v(i,j-1))
          vue(i,j)=v(i+1,j-1)+(1.d0/(1.d0+rp))*(v(i+1,j)-v(i+1,j-1))

            vu(i,j)=0.5d0*(vuw(i,j)+vue(i,j))
          rou(i,j)= (ro(i,j)+ro(i+1,j))/2.d0
          meuu(i,j)=(meu(i,j)+meu(i+1,j))/2.d0
          meuu(i,j+1)=(meu(i,j+1)+meu(i+1,j+1))/2.d0
          meuu(i,j-1)=(meu(i,j-1)+meu(i+1,j-1))/2.d0
          phiuu(i,j)=(phi(i,j)+phi(i+1,j))/2.d0

           if (u(i,j)>0.d0)then
             wu=-u(i,j)
             pu=u(i,j)
             euup=0.d0
             else
             wu=0.d0
             pu=-u(i,j)
             euup=u(i,j)
             endif

             if(vu(i,j)>0.d0)then
             sv=-vu(i,j)
             pv=vu(i,j)
             nvvp=0.d0
             else
             sv=0.d0
             pv=-vu(i,j)/rp
             nvvp=vu(i,j)/rp
             endif




          uW(i,j)= wu/dx + (1.d0/(Re0*2.d0*dx*dx*meu0))*(ro0/rou(i,j))*
     &     (meu(i+1,j)-meu(i,j))-(1.d0/(Re0*dx*dx*meu0))*(ro0/rou(i,j))*
     &      meuu(i,j)

             uE(i,j)= -(1.d0/(Re0*2.d0*dx*dx*meu0))*(ro0/rou(i,j))*
     &  (meu(i+1,j)-meu(i,j))-(1.d0/(Re0*dx*dx*meu0))*(ro0/rou(i,j))*
     &  meuu(i,j) + euup/dx

         uS(i,j)= (1.d0/(Re0*meu0))*(ro0/rou(i,j))*(meuu(i,j+1)/(rp*dy*
     & (1.d0+rp))-(1.d0-rp*rp)*meuu(i,j)/
     & (rp*dy*(1.d0+rp))-meuu(i,j-1)*rp*rp/
     & (rp*dy*(1.d0+rp)))*rp/(dy*(1.d0+rp))-
     & (1.d0/(Re0*meu0))*(ro0/rou(i,j))*
     &   meuu(i,j)*2.d0/((1.d0+rp)*dy*dy)+sv/dy

        uN(i,j)= nvvp/dy-(1.d0/(Re0*meu0))*(ro0/rou(i,j))*(meuu(i,j+1)/
     & (rp*dy*(1.d0+rp))-(1.d0-rp*rp)*meuu(i,j)/
     & (rp*dy*(1.d0+rp))-meuu(i,j-1)*
     &  rp*rp/(rp*dy*(1.d0+rp)))/(rp*(dy*(1.d0+rp)))-(1.d0/(Re0*meu0))*
     &  (ro0/rou(i,j))*meuu(i,j)*2.d0/(rp*(1.d0+rp)*dy*dy)


        uPP(i,j)=1.d0/dt+pu/dx+pv/dy+(2.d0*meuu(i,j)/(Re0*meu0*dx**2))*
     & (ro0/rou(i,j))+(2.d0*meuu(i,j)/
     & (Re0*meu0*rp*dy*dy))*(ro0/rou(i,j))+
     & (1.d0/(Re0*meu0))*(ro0/rou(i,j))*(meuu(i,j+1)/(rp*dy*(1.d0+rp))-
     &(1.d0-rp*rp)*meuu(i,j)/(rp*dy*(1.d0+rp))-meuu(i,j-1)*rp*rp/(rp*dy*
     &  (1.d0+rp)))*(1.d0-rp*rp)/(rp*dy*(1.d0+rp))

           uB(i,j)=u(i,j)/dt-(ro0/rou(i,j))*(p(i+1,j)-p(i,j))/dx+
     & Ri*(betabfu(i,j)/bitabf)*(theta(i,j)+theta(i+1,j))*0.5d0*
     & (robflu(i,j)/rou(i,j))-
     & (F)*phi0*(phiuu(i,j)-1.d0)*(rop-robflu(i,j))/rou(i,j)



            enddo
            enddo


         do i=1,N-1
         do j=2,M-1
         dterm=abs(uPP(i,j))
         nondterm=abs(uE(i,j))+abs(uW(i,j))+abs(uN(i,j))+
     &   abs(uS(i,j))
         if (dterm<nondterm)then
         flag=1
         print*,"u dia dom violated",j
         endif
         enddo
         enddo

         if(flag>0)exit





          !! coefficient for v equation
            do i=2,N-1
            do j=2,M-1                             !(M/2)-1

            if (j<(M3))then
            rv1=1.d0
            dy=dy1
            dyp=dy1
            else if(j==M3)then
            dy=dy1
            dyp=dy1
            rv1=(1.d0+(dy2/dy1))/2.d0
            else if(j==(M3+1))then
            dy=(dy1+dy2)/2.d0
            dyp=dy2
            rv1=2.d0*(dy2/dy1)/(1.d0+(dy2/dy1))

            else if(j<(M2+M3-1))then
            rv1=1.d0
            dy=dy2
            dyp=dy2
            else if(j==(M2+M3-1))then
            dy=dy2
            dyp=dy2
            rv1=(1.d0+(dy3/dy2))/2.d0
            else if(j==(M2+M3))then
            dy=(dy2+dy3)/2.d0
            dyp=dy3
            rv1=2.d0*(dy3/dy2)/(1.d0+(dy3/dy2))
            else
            dy=dy3
            rv1=1.d0
            dyp=dy3

           end if




         rov(i,j)=(ro(i,j)+ro(i,j+1))/2.d0
          meuv(i,j)=(meu(i,j)+meu(i,j+1))/2.d0
          meuv(i+1,j)=(meu(i+1,j)+meu(i+1,j+1))/2.d0
          meuv(i-1,j)=(meu(i-1,j)+meu(i-1,j+1))/2.d0
          uv(i,j)=(u(i,j)+u(i,j+1)+u(i-1,j+1)+u(i-1,j))/4.d0
           if (uv(i,j)>0.d0)then
             wuv=-uv(i,j)
             puv=uv(i,j)
             euv=0.d0
             else
             wuv=0.d0
             puv=-uv(i,j)
             euv=uv(i,j)
             endif

             if(v(i,j)>0.d0)then
             svv=-v(i,j)
             pvv=v(i,j)
             nvv=0.d0
             else
             svv=0.d0
             pvv=-v(i,j)/rv1
             nvv=v(i,j)/rv1
             endif


          vW(i,j)= wuv/dx+(1.d0/(Re0*4.d0*dx*dx*meu0))*(ro0/rov(i,j))*
     & (meuv(i+1,j)-meuv(i-1,j))-(1.d0/(Re0*dx*dx*meu0))*(ro0/rov(i,j))*
     &  meuv(i,j)

         vE(i,j)= -(1.d0/(Re0*4.d0*dx*dx*meu0))*(ro0/rov(i,j))*
     & (meuv(i+1,j)-meuv(i-1,j))-(1.d0/(Re0*dx*dx*meu0))*(ro0/rov(i,j))*
     &  meuv(i,j) + euv/dx

          vS(i,j)= (1.d0/(Re0*dy*dyp*meu0))*(ro0/rov(i,j))*(meu(i,j+1)-
     &  meu(i,j))*rv1/(1.d0+rv1)-(1.d0/(Re0*dy*dy*meu0))*(ro0/rov(i,j))*
     &    meuv(i,j)*2.d0/(1.d0+rv1) + svv/dy

          vN(i,j)= nvv/dy-(1.d0/(Re0*dy*dyp*meu0))*(ro0/rov(i,j))*
     &  (meu(i,j+1)-meu(i,j))/(rv1*(1.d0+rv1))-(1.d0/(Re0*dy*dy*meu0))*
     &  (ro0/rov(i,j))*meuv(i,j)*2.d0/(rv1*(1.d0+rv1))

        vPP(i,j)=1.d0/dt+puv/dx+pvv/dy+(2.d0*meuv(i,j)*(ro0/rov(i,j)))/
     & (meu0*Re0*dx**2)+ (1.d0*meuv(i,j)/
     &  (meu0*Re0*dy**2))*(ro0/rov(i,j))*
     &  (2.d0/(rv1))+(1.d0/(Re0*dy*dyp*meu0))*(ro0/rov(i,j))*
     & (meu(i,j+1)-meu(i,j))*(1.d0-rv1*rv1)/(rv1*(1.d0+rv1))

          vB(i,j)=v(i,j)/dt-(ro0/rov(i,j))*(p(i,j+1)-p(i,j))/dyp


         end do

         enddo
         !!image point for v at exit
         
            i=N
            do j=2,M-1                             !(M/2)-1

            if (j<(M3))then
            rv1=1.d0
            dy=dy1
            dyp=dy1
            else if(j==M3)then
            dy=dy1
            dyp=dy1
            rv1=(1.d0+(dy2/dy1))/2.d0
            else if(j==(M3+1))then
            dy=(dy1+dy2)/2.d0
            dyp=dy2
            rv1=2.d0*(dy2/dy1)/(1.d0+(dy2/dy1))

            else if(j<(M2+M3-1))then
            rv1=1.d0
            dy=dy2
            dyp=dy2
            else if(j==(M2+M3-1))then
            dy=dy2
            dyp=dy2
            rv1=(1.d0+(dy3/dy2))/2.d0
            else if(j==(M2+M3))then
            dy=(dy2+dy3)/2.d0
            dyp=dy3
            rv1=2.d0*(dy3/dy2)/(1.d0+(dy3/dy2))
            else
            dy=dy3
            rv1=1.d0
            dyp=dy3

           end if




         rov(i,j)=(ro(i,j)+ro(i,j+1))/2.d0
          meuv(i,j)=(meu(i,j)+meu(i,j+1))/2.d0
         
          meuv(i-1,j)=(meu(i-1,j)+meu(i-1,j+1))/2.d0
          uv(i,j)=(u(i,j)+u(i,j+1)+u(i-1,j+1)+u(i-1,j))/4.d0
           if (uv(i,j)>0.d0)then
             wuv=-uv(i,j)
             puv=uv(i,j)
             euv=0.d0
             else
             wuv=0.d0
             puv=-uv(i,j)
             euv=uv(i,j)
             endif

             if(v(i,j)>0.d0)then
             svv=-v(i,j)
             pvv=v(i,j)
             nvv=0.d0
             else
             svv=0.d0
             pvv=-v(i,j)/rv1
             nvv=v(i,j)/rv1
             endif


          vW(i,j)= wuv/dx+(1.d0/(Re0*4.d0*dx*dx*meu0))*(ro0/rov(i,j))*
     & (3.d0*meuv(N,j)-4.d0*meuv(N-1,j)+meuv(N-2,j))-
     & (1.d0/(Re0*dx*dx*meu0))*(ro0/rov(i,j))*meuv(i,j)

         vE(i,j)= -(1.d0/(Re0*4.d0*dx*dx*meu0))*(ro0/rov(i,j))*
     & (3.d0*meuv(N,j)-4.d0*meuv(N-1,j)+meuv(N-2,j))-
     & (1.d0/(Re0*dx*dx*meu0))*(ro0/rov(i,j))*meuv(i,j) + euv/dx

          vS(i,j)= (1.d0/(Re0*dy*dyp*meu0))*(ro0/rov(i,j))*(meu(i,j+1)-
     &  meu(i,j))*rv1/(1.d0+rv1)-(1.d0/(Re0*dy*dy*meu0))*(ro0/rov(i,j))*
     &    meuv(i,j)*2.d0/(1.d0+rv1) + svv/dy

          vN(i,j)= nvv/dy-(1.d0/(Re0*dy*dyp*meu0))*(ro0/rov(i,j))*
     &  (meu(i,j+1)-meu(i,j))/(rv1*(1.d0+rv1))-(1.d0/(Re0*dy*dy*meu0))*
     &  (ro0/rov(i,j))*meuv(i,j)*2.d0/(rv1*(1.d0+rv1))

        vPP(i,j)=1.d0/dt+puv/dx+pvv/dy+(2.d0*meuv(i,j)*(ro0/rov(i,j)))/
     & (meu0*Re0*dx**2)+ (1.d0*meuv(i,j)/
     &  (meu0*Re0*dy**2))*(ro0/rov(i,j))*
     &  (2.d0/(rv1))+(1.d0/(Re0*dy*dyp*meu0))*(ro0/rov(i,j))*
     & (meu(i,j+1)-meu(i,j))*(1.d0-rv1*rv1)/(rv1*(1.d0+rv1))

          vB(i,j)=v(i,j)/dt-(ro0/rov(i,j))*(p(i,j+1)-p(i,j))/dyp


         end do


         
          do i=2,N
         do j=1,M-1
         dterm=abs(vPP(i,j))
         nondterm=abs(vE(i,j))+abs(vW(i,j))+abs(vN(i,j))+
     &   abs(vS(i,j))
         if (dterm<nondterm)then
         flag=1
         print*,"v dia dom violated i=",i,"j=",j
         endif
         enddo
         enddo

         if(flag>0)exit

         !! calculation of temperature coefficient



          do i=2,N-1
          do j=2,M-1

          if (j<(M3+1))then
            rp=1.d0
            dy=dy1
            else if(j==(M3+1))then
            rp=dy2/dy1
            dy=dy1
            else if(j<M2+M3)then
            rp=1.d0
            dy=dy2
            else if(j==M2+M3)then
            dy=dy2
            rp=dy3/dy2
            else
            rp=1.d0
            dy=dy3
           end if

                                                        !M-1
              up(i,j)=(u(i,j)+u(i-1,j))/2.d0

              vp(i,j)=v(i,j-1)+(1.d0/(1.d0+rp))*(v(i,j)-v(i,j-1))
              if (up(i,j)>0.d0)then
             wup=-up(i,j)
             pup=up(i,j)
             eup=0.d0
             else
             wup=0.d0
             pup=-up(i,j)
             eup=up(i,j)
             endif

             if(vp(i,j)>0.d0)then
             svp=-vp(i,j)
             pvp=vp(i,j)
             nvp=0.d0
             else
             svp=0.d0
             pvp=-vp(i,j)/rp
             nvp=vp(i,j)/rp
             endif

           thW(i,j)=wup/dx+ro0*c0*(k(i+1,j)-k(i-1,j))/(Re0*Pr0*4.d0*
     &     k0*dx*dx*ro(i,j)*c(i,j))-k(i,j)*ro0*c0/(ro(i,j)*c(i,j)*
     &     Re0*Pr0*k0*dx*dx)+ro0*c0*(Db(i,j)/Db0)*
     &     (phi(i+1,j)-phi(i-1,j))/
     &   (4.d0*dx*dx*ro(i,j)*c(i,j)*Re0*Pr0*Le0)+ro0*c0*(Dtl(i,j)/Dtl0)*
     &   (theta(i+1,j)-theta(i-1,j))/(4.d0*dx*dx*ro(i,j)*c(i,j)*
     &      Re0*Pr0*Le0*Nbt0)


            thE(i,j)=-ro0*c0*(k(i+1,j)-k(i-1,j))/(Re0*Pr0*4.d0*
     &     k0*dx*dx*ro(i,j)*c(i,j))-k(i,j)*ro0*c0/(ro(i,j)*c(i,j)*
     &      Re0*Pr0*k0*dx*dx)-
     &     ro0*c0*(Db(i,j)/Db0)*(phi(i+1,j)-phi(i-1,j))/
     &   (4.d0*dx*dx*ro(i,j)*c(i,j)*Re0*Pr0*Le0)-ro0*c0*(Dtl(i,j)/Dtl0)*
     &  (theta(i+1,j)-theta(i-1,j))/(4.d0*dx*dx*ro(i,j)*c(i,j)*
     &      Re0*Pr0*Le0*Nbt0) +eup/dx

            thS(i,j)=ro0*c0*(k(i,j+1)-(1-rp*rp)*k(i,j)-rp*rp*k(i,j-1))*
     & rp/(rp*(1.d0+rp)*(1.d0+rp)*Re0*Pr0*k0*dy*dy*ro(i,j)*c(i,j))-
     & (k(i,j)*
     & ro0*c0/(ro(i,j)*c(i,j)*Re0*Pr0*k0*dy*dy))*2.d0/(1.d0+rp)+ro0*c0*
     & (Db(i,j)/Db0)*(phi(i,j+1)-(1.d0-rp*rp)*phi(i,j)-rp*rp*
     &  phi(i,j-1))*rp/(rp*(1.d0+rp)*(1.d0+rp)*
     & dy*dy*ro(i,j)*c(i,j)*Re0*Pr0*
     & Le0)+ro0*c0*(Dtl(i,j)/Dtl0)*(theta(i,j+1)-
     & (1.d0-rp*rp)*theta(i,j)-
     & rp*rp*theta(i,j-1))*rp/(rp*(1.d0+rp)*(1.d0+rp)*dy*dy*ro(i,j)*
     & c(i,j)*Re0*Pr0*Le0*Nbt0)+svp/dy

        thN(i,j)=nvp/dy-ro0*c0*(k(i,j+1)/(rp*dy*(1.d0+rp))-(1.d0-rp*rp)*
     &  k(i,j)/(rp*dy*(1.d0+rp))-rp*rp*k(i,j-1)/
     & (rp*dy*(1.d0+rp)))/(rp*(1.d0+rp)*
     &  Re0*Pr0*k0*dy*ro(i,j)*c(i,j))-(k(i,j)*ro0*c0/(ro(i,j)*c(i,j)*
     & Re0*Pr0*k0*dy*dy))*2.d0/(rp*(1.d0+rp))-ro0*c0*(Db(i,j)/Db0)*
     & (phi(i,j+1)-(1.d0-rp*rp)*phi(i,j)-
     & rp*rp*phi(i,j-1))/(rp*rp*(1.d0+rp)*
     & (1.d0+rp)*dy*dy*ro(i,j)*c(i,j)*Re0*Pr0*Le0)-ro0*c0*
     & (Dtl(i,j)/Dtl0)*
     & (theta(i,j+1)-(1.d0-rp*rp)*theta(i,j)-rp*rp*theta(i,j-1))/(rp*rp*
     & (1.d0+rp)*(1.d0+rp)*dy*dy*ro(i,j)*c(i,j)*Re0*Pr0*Le0*Nbt0)

        thP(i,j)=1.d0/dt+pup/dx+pvp/dy+ro0*c0*2.d0*k(i,j)/(Re0*k0*dx*dx*
     &    Pr0*ro(i,j)*c(i,j))+ro0*c0*2.d0*k(i,j)/(rp*Re0*k0*dy*dy*Pr0*
     &  ro(i,j)*c(i,j))+ro0*c0*(k(i,j+1)-(1.d0-rp*rp)*k(i,j)-rp*rp*
     &  k(i,j-1))*(1.d0-rp*rp)/(rp*rp*(1.d0+rp)*
     &  (1.d0+rp)*Re0*Pr0*k0*dy*dy*
     &  ro(i,j)*c(i,j))+ro0*c0*(Db(i,j)/Db0)*(phi(i,j+1)-
     &  (1.d0-rp*rp)*phi(i,j)-rp*rp*phi(i,j-1))*
     &  (1.d0-rp*rp)/(rp*rp*(1.d0+rp)*
     &  (1.d0+rp)*dy*dy*ro(i,j)*c(i,j)*Re0*Pr0*Le0)+ro0*c0*
     &  (Dtl(i,j)/Dtl0)*
     &  (theta(i,j+1)-(1.d0-rp*rp)*theta(i,j)-rp*rp*theta(i,j-1))*
     &  (1.d0-rp*rp)/(rp*rp*(1.d0+rp)*(1.d0+rp)*dy*dy*ro(i,j)*c(i,j)*
     &   Re0*Pr0*Le0*Nbt0)




          thB(i,j)=theta(i,j)/dt


          enddo


            enddo

          !!! image point for theta on exit
          i=N
          do j=2,M-1

          if (j<(M3+1))then
            rp=1.d0
            dy=dy1
            else if(j==(M3+1))then
            rp=dy2/dy1
            dy=dy1
            else if(j<M2+M3)then
            rp=1.d0
            dy=dy2
            else if(j==M2+M3)then
            dy=dy2
            rp=dy3/dy2
            else
            rp=1.d0
            dy=dy3
           end if

                                                        !M-1
              up(i,j)=(u(i,j)+u(i-1,j))/2.d0

              vp(i,j)=v(i,j-1)+(1.d0/(1.d0+rp))*(v(i,j)-v(i,j-1))
              if (up(i,j)>0.d0)then
             wup=-up(i,j)
             pup=up(i,j)
             eup=0.d0
             else
             wup=0.d0
             pup=-up(i,j)
             eup=up(i,j)
             endif

             if(vp(i,j)>0.d0)then
             svp=-vp(i,j)
             pvp=vp(i,j)
             nvp=0.d0
             else
             svp=0.d0
             pvp=-vp(i,j)/rp
             nvp=vp(i,j)/rp
             endif

           thW(i,j)=wup/dx+ro0*c0*(3.d0*k(N,j)-4.d0*k(N-1,j)+k(N-2,j))/
     &     (Re0*Pr0*4.d0*
     &     k0*dx*dx*ro(i,j)*c(i,j))-k(i,j)*ro0*c0/(ro(i,j)*c(i,j)*
     &     Re0*Pr0*k0*dx*dx)+ro0*c0*(Db(i,j)/Db0)*
     &     (3.d0*phi(N,j)-4.d0*phi(N-1,j)+phi(N-2,j))/
     &   (4.d0*dx*dx*ro(i,j)*c(i,j)*Re0*Pr0*Le0)+ro0*c0*(Dtl(i,j)/Dtl0)*
     &   (3.d0*theta(N,j)-4.d0*theta(N-1,j)+theta(N-2,j))/
     &   (4.d0*dx*dx*ro(i,j)*c(i,j)*
     &      Re0*Pr0*Le0*Nbt0)


            thE(i,j)=-ro0*c0*(3.d0*k(N,j)-4.d0*k(N-1,j)+k(N-2,j))/
     &       (Re0*Pr0*4.d0*
     &     k0*dx*dx*ro(i,j)*c(i,j))-k(i,j)*ro0*c0/(ro(i,j)*c(i,j)*
     &      Re0*Pr0*k0*dx*dx)-
     &  ro0*c0*(Db(i,j)/Db0)*(3.d0*phi(N,j)-4.d0*phi(N-1,j)+phi(N-2,j))/
     &  (4.d0*dx*dx*ro(i,j)*c(i,j)*Re0*Pr0*Le0)-ro0*c0*(Dtl(i,j)/Dtl0)*
     &  (3.d0*theta(N,j)-4.d0*theta(N-1,j)+theta(N-2,j))/
     &  (4.d0*dx*dx*ro(i,j)*c(i,j)*
     &      Re0*Pr0*Le0*Nbt0) +eup/dx

            thS(i,j)=ro0*c0*(k(i,j+1)-(1-rp*rp)*k(i,j)-rp*rp*k(i,j-1))*
     & rp/(rp*(1.d0+rp)*(1.d0+rp)*Re0*Pr0*k0*dy*dy*ro(i,j)*c(i,j))-
     & (k(i,j)*
     & ro0*c0/(ro(i,j)*c(i,j)*Re0*Pr0*k0*dy*dy))*2.d0/(1.d0+rp)+ro0*c0*
     & (Db(i,j)/Db0)*(phi(i,j+1)-(1.d0-rp*rp)*phi(i,j)-rp*rp*
     &  phi(i,j-1))*rp/(rp*(1.d0+rp)*(1.d0+rp)*
     & dy*dy*ro(i,j)*c(i,j)*Re0*Pr0*
     & Le0)+ro0*c0*(Dtl(i,j)/Dtl0)*(theta(i,j+1)-
     & (1.d0-rp*rp)*theta(i,j)-
     & rp*rp*theta(i,j-1))*rp/(rp*(1.d0+rp)*(1.d0+rp)*dy*dy*ro(i,j)*
     & c(i,j)*Re0*Pr0*Le0*Nbt0)+svp/dy

        thN(i,j)=nvp/dy-ro0*c0*(k(i,j+1)/(rp*dy*(1.d0+rp))-(1.d0-rp*rp)*
     &  k(i,j)/(rp*dy*(1.d0+rp))-rp*rp*k(i,j-1)/
     &  (rp*dy*(1.d0+rp)))/(rp*(1.d0+rp)*
     &  Re0*Pr0*k0*dy*ro(i,j)*c(i,j))-(k(i,j)*ro0*c0/(ro(i,j)*c(i,j)*
     &  Re0*Pr0*k0*dy*dy))*2.d0/(rp*(1.d0+rp))-ro0*c0*(Db(i,j)/Db0)*
     & (phi(i,j+1)-(1.d0-rp*rp)*phi(i,j)-
     & rp*rp*phi(i,j-1))/(rp*rp*(1.d0+rp)*
     & (1.d0+rp)*dy*dy*ro(i,j)*c(i,j)*Re0*Pr0*Le0)-ro0*c0*
     & (Dtl(i,j)/Dtl0)*
     & (theta(i,j+1)-(1.d0-rp*rp)*theta(i,j)-rp*rp*theta(i,j-1))/(rp*rp*
     & (1.d0+rp)*(1.d0+rp)*dy*dy*ro(i,j)*c(i,j)*Re0*Pr0*Le0*Nbt0)

        thP(i,j)=1.d0/dt+pup/dx+pvp/dy+ro0*c0*2.d0*k(i,j)/(Re0*k0*dx*dx*
     &    Pr0*ro(i,j)*c(i,j))+ro0*c0*2.d0*k(i,j)/(rp*Re0*k0*dy*dy*Pr0*
     &  ro(i,j)*c(i,j))+ro0*c0*(k(i,j+1)-(1.d0-rp*rp)*k(i,j)-rp*rp*
     & k(i,j-1))*(1.d0-rp*rp)/(rp*rp*(1.d0+rp)*
     & (1.d0+rp)*Re0*Pr0*k0*dy*dy*
     & ro(i,j)*c(i,j))+ro0*c0*(Db(i,j)/Db0)*(phi(i,j+1)-
     & (1.d0-rp*rp)*phi(i,j)-rp*rp*phi(i,j-1))*
     &  (1.d0-rp*rp)/(rp*rp*(1.d0+rp)*
     &(1.d0+rp)*dy*dy*ro(i,j)*c(i,j)*Re0*Pr0*Le0)+ro0*c0*
     &  (Dtl(i,j)/Dtl0)*
     &  (theta(i,j+1)-(1.d0-rp*rp)*theta(i,j)-rp*rp*theta(i,j-1))*
     & (1.d0-rp*rp)/(rp*rp*(1.d0+rp)*(1.d0+rp)*dy*dy*ro(i,j)*c(i,j)*
     &  Re0*Pr0*Le0*Nbt0)




          thB(i,j)=theta(i,j)/dt


          enddo







          !! coefficient calculation completed

          do i=2,N
         do j=2,M-1
         dterm=abs(thP(i,j))
         nondterm=abs(thE(i,j))+abs(thW(i,j))+abs(thN(i,j))+
     &   abs(thS(i,j))
         if (dterm<nondterm)then
         flag=1
         print*,"theta dia dom violated"
         endif
         enddo
         enddo

         if(flag>0)exit

          do i=1,N
          do j=1,M
          phiold(i,j)=phi(i,j)
          enddo
          enddo


         
           phipe=30.d0

           do while (phipe>(0.0000001d0*dt))
            
            phiu=0.d0
            

           
          do i=1,N+1
          do j=0,M
           phieu(i,j)=phi(i,j)
          enddo
          enddo

          do i=2,N
           do j=1,M-1
          phi(i,j)=0.7d0*((phiB(i,j)-(phiW(i,j)*phi(i-1,j)+phiE(i,j)*
     &   phi(i+1,j)+phiS(i,j)*phi(i,j-1)+phiN(i,j)*
     &   phi(i,j+1)))/phiP(i,j))+0.3d0*phi(i,j)

           end do
          end do

            do i=2,N
         phi(i,0)=phi(i,2)+K1(i)
         phi(i,M)=1.d0

            enddo

            do j=1,M
             phi(1,j)=1.d0
            enddo
            do j=1,M-1
             phi(N+1,j)=phi(N-1,j)
            end do


           do i=2,N
           do j=0,M-1
           phir(i,j)=ABS(phieu(i,j)-phi(i,j))
           if(abs(phir(i,j))>phiu)then
           phiu=abs(phir(i,j))
           endif
           enddo
           end do
            phipe=phiu
          
           end do !! phi equation closed
          

           do i=0,N
           do j=1,M
           uuoold(i,j)=u(i,j)
           enddo
           enddo
             !!! u momentum equation starts
                !! u equation starts
           upe=30.d0

           do while (upe>(0.000001d0*dt))
          
            uiu=0.d0


           do i=0,N
           do j=1,M
           eu(i,j)=u(i,j)
           end do
           enddo

           do i=1,N-1
           do j=2,M-1

          u(i,j)=0.7d0*(1/uPP(i,j))*(uB(i,j)-(uE(i,j)*u(i+1,j)+
     &    uW(i,j)*u(i-1,j)+uS(i,j)*u(i,j-1)+uN(i,j)*u(i,j+1)))+
     &    0.3d0*u(i,j)




           end do
          end do

          do i=0,N
           u(i,1)=1.d0
           u(i,M)=0.d0
          enddo

          do j=2,M-1
            u(0,j)=0.d0
            u(N,j)=u(N-1,j)
          enddo

           do i=1,N-1
           do j=1,M-1
           ur(i,j)=ABS(eu(i,j)-u(i,j))
           if(abs(ur(i,j))>uiu)then
           uiu=abs(ur(i,j))
           endif
           enddo
           end do
            upe=uiu
           
           !print*,"uiu",uiu
           end do !! u equation closed

          

             do i=1,N
           do j=0,M
           vold(i,j)=v(i,j)
           enddo
           enddo

               !! v equation starts
           vpe=30.d0

           do while (vpe>(0.0000001d0*dt))
           
            viu=0.d0


          do i=1,N+1
           do j=0,M
           ev(i,j)=v(i,j)
           enddo
           enddo

           do i=2,N
           do j=1,M-1

          v(i,j)=0.7d0*(vB(i,j)-(vW(i,j)*v(i-1,j)+vE(i,j)*v(i+1,j)+
     &    vS(i,j)*v(i,j-1)+vN(i,j)*v(i,j+1)))/vPP(i,j)+0.3d0*v(i,j)

           end do
          end do

          do i=1,N
            v(i,0)=0.d0
          enddo
          do i=2,N
            v(i,M)=v(i,M-1)
          enddo

          do j=1,M-1
            v(1,j)=0.d0
           v(N+1,j)=v(N-1,j)
          enddo





           do i=2,N
           do j=1,M-1
           vr(i,j)=ABS(ev(i,j)-v(i,j))
           if(abs(vr(i,j))>viu)then
           viu=abs(vr(i,j))
           endif
           enddo
           end do
            vpe=viu
           end do !! v equation closed


           do i=1,N
           do j=1,M
           thold(i,j)=theta(i,j)
           enddo
           enddo

                !! T equation starts
           thpe=30.d0

           do while (thpe>(0.0000001d0*dt))
            thiu=0.d0


          do i=1,N+1
           do j=1,M
           eth(i,j)=theta(i,j)

           enddo
           enddo

            do i=2,N
           do j=2,M-1


          theta(i,j)=0.7*(thB(i,j)-(thW(i,j)*theta(i-1,j)+
     &     thE(i,j)*theta(i+1,j)+thS(i,j)*theta(i,j-1)+
     &     thN(i,j)*theta(i,j+1)))/thP(i,j)+0.3*theta(i,j)

           end do
          end do

           do i=1,N
        theta(i,1)=1.d0
        theta(i,M)=0.d0
           enddo

           do j=2,M-1
            theta(1,j)=0.d0
           theta(N+1,j)=theta(N-1,j)
           enddo




           do i=2,N
           do j=1,M-1
           thr(i,j)=ABS(eth(i,j)-theta(i,j))
           if(abs(thr(i,j))>thiu)then
           thiu=abs(thr(i,j))
           endif
           enddo
           end do
            thpe=thiu
          end do !! T equation closed

         !!!!!!!!! pressure correction equation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !==============================================================================




            do i=1,N
            do j=1,M
            thetaee(i,j)=theta(i,j)-thold(i,j)
            phieee(i,j)=phi(i,j)-phiold(i,j)
            enddo
            enddo

        !!! calculating coefficient for pressure correction equation

           epp(1,1)=0.d0
           epp(1,M)=0.d0
         do i=2,N-1
          epp(i,1)=(-v(i,2)+9.d0*v(i,1)-8.d0*v(i,0))/(3.d0*dy1)
          !epp(i,M)=(v(i,M)-v(i,M-1))/(dy3)
         enddo

         do j=2,M-1
           epp(1,j)=(-u(2,j)+9.d0*u(1,j)-8.d0*u(0,j))/(3.d0*dx)
         enddo

          do i=2,N-1
           do j=2,M3                                                     !M-1

            dy=dy1
          epp(i,j)=(u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i,j-1))/dy
          enddo

          j=M3+1
          rp=dy2/dy1
          dy=dy1
           vp(i,j)=v(i,j-1)+(1/(1+rp))*(v(i,j)-v(i,j-1))
            vp(i,j+1)=0.5d0*(v(i,j)+v(i,j+1))
            vp(i,j-1)=0.5d0*(v(i,j-2)+v(i,j-1))
          epp(i,j)=(u(i,j)-u(i-1,j))/dx+(vp(i,j+1)/(rp*(1.d0+rp))-
     &   (1.d0-rp*rp)*vp(i,j)/(rp*(1.d0+rp))-vp(i,j-1)*rp/(1.d0+rp))/dy



          do j=M3+2,M2+M3-1

              dy=dy2
          epp(i,j)=(u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i,j-1))/dy

          enddo
            j=M2+M3
             rp=dy3/dy2
          dy=dy2
           vp(i,j)=v(i,j-1)+(1.d0/(1.d0+rp))*(v(i,j)-v(i,j-1))
            vp(i,j+1)=0.5d0*(v(i,j)+v(i,j+1))
            vp(i,j-1)=0.5d0*(v(i,j-2)+v(i,j-1))
          epp(i,j)=(u(i,j)-u(i-1,j))/dx+(vp(i,j+1)/(rp*(1.d0+rp))-
     &   (1.d0-rp*rp)*vp(i,j)/(rp*(1.d0+rp))-vp(i,j-1)*rp/(1.d0+rp))/dy

           do j=M2+M3+1,M-1
              dy=dy3
          epp(i,j)=(u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i,j-1))/dy

          enddo
          enddo

            do i=1,N-1
           do j=1,M-1
             !============================================================================================
            if (j<(M3+1))then
            rp=1.d0
            dy=dy1
            else if(j==(M3+1))then
            rp=dy2/dy1
            dy=dy1
            else if(j<M2+M3)then
            rp=1.d0
            dy=dy2
            else if(j==M2+M3)then
            dy=dy2
            rp=dy3/dy2
            else
            rp=1.d0
            dy=dy3
           end if
          if(i==1)then
           peu(1,j)=(epp(1,j))/A  +
     &  Ri*(betabf(1,j)/bitabf)*(robfl(1,j)/ro(1,j))*
     & (4.d0*thetaee(2,j)-thetaee(3,j)-3.d0*thetaee(1,j))/(2.d0*dx)-
     & (F)*((rop-robfl(1,j))/ro(1,j))*phi0*
     & (4.d0*phieee(2,j)-phieee(3,j)-3.d0*phieee(1,j))/(2.d0*dx)
          pepx(1,j)=1.d0/(dx*dx)
          pepys(1,j)=2.d0/((1.d0+rp)*(dy*dy))
          pepyn(1,j)= 2.d0/(rp*(1.d0+rp)*dy*dy)
          pep(1,j)=-2.d0/(dx*dx)-2.d0/(rp*dy*dy)

          
          else if(i>1.and.i<N)then
                                                         !M-1
           peu(i,j)=(epp(i,j))/A   +
     &  Ri*(betabf(i,j)/bitabf)*(robfl(i,j)/ro(i,j))*
     &  (thetaee(i+1,j)-thetaee(i-1,j))/(2.d0*dx)-(F)*
     & ((rop-robfl(i,j))/ro(i,j))*phi0*(phieee(i+1,j)-
     &  phieee(i-1,j))/(2.d0*dx)
          pepx(i,j)=1.d0/(dx*dx)
          pepys(i,j)=2.d0/((1.d0+rp)*(dy*dy))
          pepyn(i,j)= 2.d0/(rp*(1.d0+rp)*dy*dy)
          pep(i,j)=-2.d0/(dx*dx)-2.d0/(rp*dy*dy)



         endif

          enddo
          enddo


            err=30.d0
            t=0

         do while (err>(0.0000001*dt))
          !do t=1,100



          ! pressure correction equation



          do i=0,N
           do j=0,M
            pee(i,j)=pe(i,j)
            enddo
           enddo

          do i=1,N-1
           do j=1,M-1
         pe(i,j)=0.7d0*(peu(i,j)-(pepx(i,j)*pe(i-1,j)+
     &    pepx(i,j)*pe(i+1,j)+
     &  pepys(i,j)*pe(i,j-1)+pepyn(i,j)*pe(i,j+1)))/pep(i,j)+
     &    0.3d0*pe(i,j)
           end do
          end do

          do i=1,N-1
          pe(i,0)=pe(i,2)
          pe(i,M)=0.d0
          enddo

          do j=1,M-1
          pe(0,j)=pe(2,j)
          enddo
		  
          do j=1,M
          pe(N,j)=0.d0
          enddo


          e=0.d0

          do i=1,N-1
           do j=1,M-1
           er(i,j)=ABS(pee(i,j)-pe(i,j))
           if(abs(er(i,j))>e)then
           e=abs(er(i,j))
           endif
          
           enddo
           enddo

           err=e

          !print*, "pressure error",e
         t=t+1

          end do


          !print*,'pe',pe(10,10)
          do i=1,N-1
           do j=1,M-1
                p(i,j)=p(i,j)+1.0*pe(i,j)
           enddo
          enddo

          do i=1,N-1
          do j=2,M-1

                u(i,j)=u(i,j)-A*(pe(i+1,j)-pe(i,j))/dx  +
     &   Ri*A*(betabfu(i,j)/bitabf)*(robflu(i,j)/rou(i,j))*
     &    (thetaee(i+1,j)+thetaee(i,j))/(2.d0)-
     &   (F)*A*((rop-robflu(i,j))/rou(i,j))*phi0*
     &    (phieee(i+1,j)+phieee(i,j))/(2.d0)
          enddo
          enddo

          do i=2,N-1
          do j=1,M-1
          if(j<(M3+1))then
          dy=dy1
          else if(j<(M2+M3))then
          dy=dy2
          else
          dy=dy3
          endif
                v(i,j)=v(i,j)-A*(pe(i,j+1)-pe(i,j))/dy

           end do
          end do
           !! checking steadiness
           ueebig=0.d0
          us_sum=0.d0
          us_num=0
           do i=1,N-1
            do j=2,M-1
              if(abs(u(i,j)-uuoold(i,j))>ueebig)then
              ueebig= abs(u(i,j)-uuoold(i,j))
              endif
              us_sum=us_sum+(abs(u(i,j)-uuoold(i,j)))**2
              us_num=us_num+1
            enddo
            enddo
            u_rmse=sqrt(us_sum/us_num)
            !if(ueebig<(0.0005*dt))exit
             veebig=0.d0
          vs_sum=0.d0
          vs_num=0
           do i=2,N
            do j=1,M-1
              if(abs(v(i,j)-vold(i,j))>veebig)then
              veebig= abs(v(i,j)-vold(i,j))
              endif
              vs_sum=vs_sum+(abs(v(i,j)-vold(i,j)))**2
              vs_num=vs_num+1
            enddo
            enddo
            v_rmse=sqrt(vs_sum/vs_num)

            phibig=0.d0
            phis_sum=0.d0
            phis_num=0
           do i=2,N
            do j=1,M-1
              if(abs(phi(i,j)-phiold(i,j))>phibig)then
              phibig= abs(phi(i,j)-phiold(i,j))
              endif
             phis_sum=phis_sum+(abs(phi(i,j)-phiold(i,j)))**2
              phis_num=phis_num+1
            enddo
            enddo
            phi_rmse=sqrt(phis_sum/phis_num)

             thbig=0.d0
             thetas_sum=0.d0
             thetas_num=0
           do i=2,N
            do j=2,M-1
              if(abs(theta(i,j)-thold(i,j))>thbig)then
              thbig= abs(theta(i,j)-thold(i,j))
              endif
              thetas_sum=thetas_sum+(abs(theta(i,j)-thold(i,j)))**2
              thetas_num=thetas_num+1
            enddo
            enddo
            theta_rmse=sqrt(thetas_sum/thetas_num)

            tolerance=u_rmse

            if(v_rmse>tolerance)then
           tolerance=v_rmse
           endif

           if(theta_rmse>tolerance)then
           tolerance=theta_rmse
           endif
           if(phi_rmse>tolerance)then
           tolerance=phi_rmse
           endif



            tolerance1=ueebig

           if(veebig>tolerance1)then
           tolerance1=veebig
           endif


           if(thbig>tolerance1)then
           tolerance1=thbig
           endif
           if(phibig>tolerance1)then
           tolerance1=phibig
           endif


            do j=1,M
             u(N,j)=u(N-1,j)
            enddo
            
            do i=1,N
             v(i,M)=v(i,M-1)
            enddo




             ki=ki+1
           print*,"is",ki

           !print*,"k1",K1(5),K1(2),K1(10)
           !print*,"pressure iteration",t        !!, "theta_rmse",theta_rmse
         
           !print*,"u_rmse",u_rmse,"phi_rmse",phi_rmse
           !print*,"dt", dt, "uebig",ueebig
           !print*,"Grnf",Grnf,"Rinf",Rinf,"Ri",Ri,"Grbf",Grbf,"Re0",Re0
           !print*,"Rebf",Reb,"U0",U0
           !print*,"thetabig",thbig
           !print*,'dy',dy
            !print*,"Db",Db(5,5),"Dt",Dtl(5,5),"k",k(5,5),"c",c(5,5)
            !print*,"k0",k0,"Db",Db0,"Dt",Dtl0,"Nbt",Nbt0,"phi",phi(5,5)
           ! print*,"T",T0+(Ts-T0)*theta(10,2)
             !print*,"t",tolerance,"t1",tolerance1
            enddo !! iteraion loop closed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(flag>0)then
           print*,"diagonal dominance violated",dx
            else
            exit
            endif
              checkn=checkn+1
            enddo
           dy=dy1
          do i=1,N
       Nul(i)=-dx*(i-1)*(-25.d0*theta(i,1)+48.d0*theta(i,2)-
     &  36.d0*theta(i,3)+
     &   16.d0*theta(i,4)-3.d0*theta(i,5))/(12.d0*dy)
       GrL(i)=betabf(i,1)*g*(Ts-T0)*((dx*(i-1))**3)*((ro(i,1))**2)/
     &  ((meu(i,1))**2)
          ReL(i)=U0*dx*(i-1)*L*ro(i,1)/meu(i,1)
          RiL(i)=GrL(i)/((ReL(i))**2)
         enddo
         sum=Nul(1)+Nul(N)
         do i=2,N-1
         sum=sum+2*Nul(i)
         enddo
         Nuavg=dx*sum/2.d0

         print*,"Nuavg",Nuavg
         print*,"tol_rmse",tolerance,"tol_max",tolerance1


         !!! calculating u,v value at scalar grid points

        do j=1,M
          us_grid(1,j)=u(0,j)
          do i=2,N
          us_grid(i,j)=(u(i,j)+u(i-1,j))/2.d0
          enddo
        enddo

        do i=1,N
         vs_grid(i,1)=v(i,0)
         do j=2,M3
          vs_grid(i,j)=(v(i,j)+v(i,j-1))/2.d0
         enddo
         j=M3+1
         rp=dy2/dy1
           vs_grid(i,j)=v(i,j-1)+(1/(1+rp))*(v(i,j)-v(i,j-1))
         do j=M3+2,M2+M3-1
         vs_grid(i,j)=(v(i,j)+v(i,j-1))/2.d0
         enddo
         j=M2+M3
         rp=dy3/dy2
         vs_grid(i,j)=v(i,j-1)+(1.d0/(1.d0+rp))*(v(i,j)-v(i,j-1))
         do j=M2+M3+1,M
          vs_grid(i,j)=(v(i,j)+v(i,j-1))/2.d0

         enddo

        enddo

        !!! creation of y grid points
        ys_grid(1)=0.d0
        do j=2,M3+1
         ys_grid(j)=ys_grid(j-1)+dy1
        enddo
        do j=M3+2,M2+M3
         ys_grid(j)=ys_grid(j-1)+dy2
        enddo
        do j=M2+M3+1,M
         ys_grid(j)=ys_grid(j-1)+dy3
        enddo

         do i=1,N
         si(i,1)=0
         enddo
         do j=1,M
         si(1,j)=0
         enddo

         do i=2,N
         do j=2,M
         si(i,j)=si(i-1,j)-0.5d0*(vs_grid(i,j)+vs_grid(i-1,j))*dx
         
         enddo
         enddo

        open(unit=5,file="u.txt")
         do i=1,N
         write(5,10)(us_grid(i,j),j=1,M)
         enddo
  10     format(1x,5000f16.6)

         open(unit=6,file="v.txt")
         do i=1,N
         write(6,11)(vs_grid(i,j),j=1,M)
         end do
  11     format(1x,5000f16.6)

         open(unit=7,file="theta.txt")
         do i=1,N
         write(7,12)(theta(i,j),j=1,M)
         enddo
  12     format(1x,5000f20.10)

         open(unit=8,file="p.txt")
         do i=1,N
         write(8,13)(p(i,j),j=1,M)
         enddo
  13     format(1x,5000f20.6)

         open(unit=9,file="pe.txt")
         do i=1,N
         write(9,14)(pe(i,j),j=1,M)
         enddo
  14     format(1x,5000f16.10)

         open(unit=10,file="si.txt")
         do i=1,N
         write(10,15)(si(i,j),j=1,M)
         enddo
  15     format(1x,5000f20.10)

         open(unit=11,file="phi.txt")
         do i=1,N
         write(11,16)(phi(i,j),j=1,M)
         enddo
  16     format(1x,5000f20.10)


         open(unit=12,file="Nuavg.txt")
         write(12,17)(Nuavg)
   17    format(1x,100f20.10)


          open(unit=13,file="Ri.txt")
         write(13,18)(Ri)
   18    format(1x,100f20.10)

          open(unit=14,file="Re.txt")
         write(14,19)(Re0)
   19    format(1x,100f20.10)

         open(unit=15,file="Pr.txt")
         write(15,20)(Pr0)
   20    format(1x,100f20.10)

          open(unit=16,file="Sc.txt")
         write(16,21)(Sc0)
   21    format(1x,100f20.10)

         open(unit=17,file="Nbt.txt")
         write(17,22)(Nbt0)
   22    format(1x,100f20.10)

         open(unit=18,file="F.txt")
         write(18,23)(F)
   23    format(1x,100f20.10)

         open(unit=19,file="LeDbDtsc.txt")
         write(19,24)Le0,Db0,Dtl0*T0,Sc0
   24    format(1x,100f40.20,1x,100f30.20,1x,100f30.20,1x,100f30.20)

         open(unit=20,file="y.txt")
         do j=1,M
         write(20,25)(ys_grid(j))
         enddo
   25    format(1x,5000f20.10)

         open(unit=21,file="GrbfRebU0.txt")
         write(21,26)Grbf,Reb,U0
   26    format(1x,100f30.10)

         open(unit=22,file="Grnf.txt")
         write(22,27)(Grnf)
   27    format(1x,100f30.5)

         open(unit=23,file="RinfRibf.txt")
         write(23,28)Rinf,Ribf
   28    format(1x,f30.10,1x,f30.10)

         open(unit=24,file="knf.txt")
         do i=1,N
         write(24,29)(k(i,j),j=1,M)
         enddo
   29    format(1x,5000f20.10)

          open(unit=25,file="RiL.txt")
         do i=1,N
         write(25,30)(RiL(i))
         enddo
   30    format(1x,5000f20.10)

         open(unit=26,file="ReLNuLFopGrLFrp.txt")
         do i=1,N
         write(26,31) ReL(i),NuL(i),NuL(i)/sqrt(ReL(i)),
     &    GrL(i),NuL(i)*((Grl(i))**(-0.25d0))
         enddo
   31    format(1x,f30.10,1x,f30.10,1x,f30.10,1x,f30.10,1x,f30.10)

          open(unit=27,file="meu.txt")
         do i=1,N
         write(27,32)(meu(i,j),j=1,M)
         enddo
  32     format(1x,5000f16.10)



         end program smple_cmp

