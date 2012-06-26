      PROGRAM main
      PARAMETER(pi=3.14159,nzx=300,NBX=22,nkx=35,NC=4,NV=8)
      REAL L0,L1,L2,L3,sz(nzx),sk(nkx),
     &     bandc(NC,nkx,nkx),bandv(NV,nkx,nkx)
      COMPLEX pcv(NC,NV,3,nkx,nkx),pcc(NC,NC,3,nkx,nkx),
     &        pvv(NV,NV,3,nkx,nkx),occ(NC,NC,nkx,nkx),
     &        ovv(NV,NV,nkx,nkx),ocv(NC,NV,nkx,nkx),ovc(NV,NC,nkx,nkx)
      EXTERNAL H110BE,H110,P110
      COMMON/block01/hb,xm0,beta
      COMMON/block02/e0,e1,d0,d1,dm,p0,p1,q,
     &               gLc,gL1,gL2,gL3,gc,g1,g2,g3,ck
      COMMON/block03/e0l,e1l,d0l,d1l,dml,p0l,p1l,ql,
     &               gLcl,gL1l,gL2l,gL3l,gcl,g1l,g2l,g3l,ckl
      COMMON/block04/e0r,e1r,d0r,d1r,dmr,p0r,p1r,qr,
     &               gLcr,gL1r,gL2r,gL3r,gcr,g1r,g2r,g3r,ckr
      COMMON/block05/L0,L1,L2,L3
      COMMON/block06/nz,sz,dz
      COMMON/block07/nk,sk,dk
      COMMON/block08/bandc,bandv
      COMMON/block09/pcc,pvv,pcv
      COMMON/block10/phi,eE0,tdu,eEz
      COMMON/block11/occ,ovv,ocv,ovc

      hb=.658229                !Planck constant (eV.fs)
      xm0=9.11*10/(1.783*9.)    !free electron mass (eV.(fs/nm)^2)
      beta=.5*(hb**2)/xm0       !(eV.nm^2)
      tem=0.1                   !K

      CALL param_AlGaAs(0.,tem,e0,e1,d0,d1,dm,p0,p1,q,
     &                  gLc,gL1,gL2,gL3,gc,g1,g2,g3,ck)
      CALL param_AlGaAs(0.35,tem,e0l,e1l,d0l,d1l,dml,p0l,p1l,ql,
     &                  gLcl,gL1l,gL2l,gL3l,gcl,g1l,g2l,g3l,ckl)
      CALL param_AlGaAs(0.35,tem,e0r,e1r,d0r,d1r,dmr,p0r,p1r,qr,
     &                  gLcr,gL1r,gL2r,gL3r,gcr,g1r,g2r,g3r,ckr)

      L1=12  !nm
      L2=12  !nm
      L3=12  !nm
      L0=L1+L2+L3  !nm
      nz=nzx
      dz=L0/(nz-1)
      DO i=1,nz
         sz(i)=(i-1)*dz-L0/2
      ENDDO

      nk=25
      skmax=.6   !nm^-1
      dk=2*skmax/(nk-1)
      DO i=1,nk
         sk(i)=(i-(nk+1)/2)*dk
      ENDDO

      CALL bandstr(H110BE,H110,P110)
      eg0=bandc(1,(nk+1)/2,(nk+1)/2)-bandv(1,(nk+1)/2,(nk+1)/2)
      WRITE(*,*)'Band gap: ',eg0

      eE0=.005      !Pulse amplitude (eV/nm)
      phi=0.        !relative phase of EMW (rad.)
      tdu=50.       !Gaussian standard deviation (fs)
      ece=0.05      !photon excess energy (eV)
      eEz=0.006		!		(ev/nm)
      OPEN(UNIT=1,FILE='cohphi.txt')
c      DO WHILE (phi.le.(2*pi+.00001))
         CALL photocurrent(eg0,ece,cohccint_x,cohccint_y,
     &        cohvvint_x,cohvvint_y,cohcvint_x,cohcvint_y)
         WRITE(1,'(22f23.16)')phi/pi,cohccint_x,cohccint_y,
     &         cohvvint_x,cohvvint_y,cohcvint_x,cohcvint_y
         WRITE(*,*)'phase=',phi,cohvvint_y,cohcvint_y
c         phi=phi+.1*pi
c      ENDDO
      CLOSE(1)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION eE(t,om,tdu,eE0,phi)   !(eV/nm)
      eE=eE0*exp(-.5*(t/tdu)**2)*cos(om*t+phi)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION eA(t,om,tdu,eE0,phi)
      PARAMETER (pi=3.14159)
      dt=0.01*(2*pi/om)
      t1=-5*tdu
      sum1=0.
      DO WHILE (t1.le.t)
         sum1=sum1+dt*eE(t1,om,tdu,eE0,phi)
         t1=t1+dt
      ENDDO
      eA=-sum1
      END
C!the vecto theo oz khi co dien truong ngoai theo phong oz
      FUNCTION eAz(t,om,tdu,eEz)   !(eV*fs/nm)
	      sumz=0.0
	      dt=0.01*(2*pi/om)
	      t1=-5*tdu
      DO WHILE (t1.le.t)
         sumz=sumz+dt*eEz
         t1=t1+dt
      ENDDO
      eAz=-sumz
      END
!frequence rabi
      SUBROUTINE rabifreq(t,om)
      PARAMETER(nkx=35,NC=4,NV=8)
      REAL sk(nkx)
      COMPLEX pcv(NC,NV,3,nkx,nkx),pcc(NC,NC,3,nkx,nkx),
     &        pvv(NV,NV,3,nkx,nkx),occ(NC,NC,nkx,nkx),
     &        ovv(NV,NV,nkx,nkx),ocv(NC,NV,nkx,nkx),ovc(NV,NC,nkx,nkx)
      COMMON/block01/hb,xm0,beta
      COMMON/block07/nk,sk,dk
      COMMON/block09/pcc,pvv,pcv
      COMMON/block10/phi,eE0,tdu,eEz
      COMMON/block11/occ,ovv,ocv,ovc

      Ax=eA(t,om,tdu,eE0,0.)/xm0               !nm/fs
      Ay=eA(t,om,tdu,eE0,phi)/xm0              !nm/fs
      Az=eAz(t,om,tdu,eEz)/xm0

      DO i=1,nk
      DO j=1,nk
         DO L1=1,NC
         DO L2=1,NC
            occ(L1,L2,i,j)=Ax*pcc(L1,L2,1,i,j)
     &                    +Ay*pcc(L1,L2,2,i,j)
     &                    +Az*pcc(L1,L2,3,i,j)
         ENDDO
         ENDDO
         DO L1=1,NV
         DO L2=1,NV
            ovv(L1,L2,i,j)=Ax*pvv(L1,L2,1,i,j)
     &                    +Ay*pvv(L1,L2,2,i,j)
     &                    +Az*pvv(L1,L2,3,i,j)
         ENDDO
         ENDDO
         DO Lc=1,NC
         DO Lv=1,NV
            ocv(Lc,Lv,i,j)=Ax*pcv(Lc,Lv,1,i,j)
     &                    +Ay*pcv(Lc,Lv,2,i,j)
     &                    +Az*pcv(Lc,Lv,3,i,j)
            ovc(Lv,Lc,i,j)=conjg(ocv(Lc,Lv,i,j))
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      END
! injection current+population charge current

      SUBROUTINE photocurrent(eg0,ece,cohccint_x,cohccint_y,
     &           cohvvint_x,cohvvint_y,cohcvint_x,cohcvint_y)
      PARAMETER (nkx=35,NC=4,NV=8,pi=3.14159,ntx=80000,nfx=10000)
      REAL sk(nkx),ye(NC,nkx,nkx),yh(NV,nkx,nkx),
     &     tme(ntx),cure_x(ntx),cure_y(ntx),curh_x(ntx),curh_y(ntx),
     &     cohcc_x(ntx),cohcc_y(ntx),cohvv_x(ntx),cohvv_y(ntx),
     &     cohcv_x(ntx),cohcv_y(ntx),
     &		omf(nfx),tmef(ntx),	! cac so hang khi bien doi fourier
     &		cohccif_x(ntx),cohccif_y(ntx),cohvvif_x(ntx),cohvvif_y(ntx),
     &		cohcvif_x(ntx),cohcvif_y(ntx)

      COMPLEX yc(NC,NC,nkx,nkx),yv(NV,NV,nkx,nkx),yp(NC,NV,nkx,nkx),
     &        pcc(NC,NC,3,nkx,nkx),pvv(NV,NV,3,nkx,nkx),
     &        pcv(NC,NV,3,nkx,nkx),sum3,sum4,sum5,
     &		cohccf_x(nfx),cohccf_y(nfx),cohvvf_x(nfx),cohvvf_y(nfx),
     &		cohcvf_x(nfx),cohcvf_y(nfx)
      COMMON/block01/hb,xm0,beta
      COMMON/block07/nk,sk,dk
      COMMON/block09/pcc,pvv,pcv
      tmin=-300.  !fs
      tmax=1000.   !fs
      t=tmin
      dt=.025     !fs
      nt=(tmax-tmin)/dt
      IF (nt.gt.ntx) STOP 'too many t-points: nt > ntx'
      om=(eg0+ece)/hb     !Photon frequency (fs-1)
! reset n,p
      DO i=1,nk
      DO j=1,nk
         DO L1=1,NC
         DO L2=1,NC
            IF (L1.eq.L2) THEN
               ye(L1,i,j)=0.
               yc(L1,L1,i,j)=0.
            ELSE
               yc(L1,L2,i,j)=(0.,0.)
            ENDIF
         ENDDO
         ENDDO
         DO L1=1,NV
         DO L2=1,NV
            IF (L1.eq.L2) THEN
               yh(L1,i,j)=0.
               yv(L1,L1,i,j)=1.
            ELSE
               yv(L1,L2,i,j)=(0.,0.)
            ENDIF
         ENDDO
         ENDDO
         DO Lc=1,NC
         DO Lv=1,NV
            yp(Lc,Lv,i,j)=(0.,0.)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
! SBE solving
      OPEN(UNIT=3,FILE='evol.txt')
      tmp1=0.
      tmp2=0.
      tmp3=0.
      tmp4=0.
      tmp5=0.
      tmp6=0.
      DO it=1,nt
         CALL RK2SBE(t,dt,om,ye,yh,yc,yv,yp,ye,yh,yc,yv,yp)
         CALL density(ye,yh,eden,hden)
         DO id=1,2
            sum1=0.
            sum2=0.
            sum3=(0.,0.)
            sum4=(0.,0.)
            sum5=(0.,0.)
            DO i=1,nk
            DO j=1,nk
               DO L=1,NC
                  sum1=sum1+(dk**2)*pcc(L,L,id,i,j)*ye(L,i,j)
               ENDDO
               DO L=1,NV
                  sum2=sum2-(dk**2)*pvv(L,L,id,i,j)*yh(L,i,j)
               ENDDO
               DO Lc=1,NC
               DO Lv=1,NV
                  sum3=sum3+(dk**2)*pcv(Lc,Lv,id,i,j)*yp(Lc,Lv,i,j)
               ENDDO
               ENDDO
               DO L1=1,NC-1
               DO L2=L1+1,NC
                  sum4=sum4+(dk**2)*pcc(L1,L2,id,i,j)*yc(L1,L2,i,j)
               ENDDO
               ENDDO
               DO L1=1,NV-1
               DO L2=L1+1,NV
                  sum5=sum5+(dk**2)*pvv(L1,L2,id,i,j)*yv(L1,L2,i,j)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
            IF (id.eq.1) THEN
               cure_x(it)=sum1/((2*pi)**2)
               curh_x(it)=sum2/((2*pi)**2)
               cohcv_x(it)=2*real(sum3)/((2*pi)**2)
               cohcc_x(it)=2*real(sum4)/((2*pi)**2)
               cohvv_x(it)=2*real(sum5)/((2*pi)**2)
            ELSE
               cure_y(it)=sum1/((2*pi)**2)
               curh_y(it)=sum2/((2*pi)**2)
               cohcv_y(it)=2*real(sum3)/((2*pi)**2)
               cohcc_y(it)=2*real(sum4)/((2*pi)**2)
               cohvv_y(it)=2*real(sum5)/((2*pi)**2)
            ENDIF
         ENDDO
         tme(it)=t
         WRITE(3,'(22f23.16)')t,eden,hden,cure_x(it),cure_y(it),
     &         curh_x(it),curh_y(it),cohcc_x(it),cohcc_y(it),
     &         cohvv_x(it),cohvv_y(it),cohcv_x(it),cohcv_y(it)
         tmp1=tmp1+dt*cohcc_x(it)
         tmp2=tmp2+dt*cohcc_y(it)
         tmp3=tmp3+dt*cohvv_x(it)
         tmp4=tmp4+dt*cohvv_y(it)
         tmp5=tmp5+dt*cohcv_x(it)
         tmp6=tmp6+dt*cohcv_y(it)
         t=t+dt
      ENDDO
      CLOSE(3)
      cohccint_x=tmp1
      cohccint_y=tmp2
      cohvvint_x=tmp3
      cohvvint_y=tmp4
      cohcvint_x=tmp5
      cohcvint_y=tmp6
c      END
!cccccccccccccccccccccc
      omfmin=0.075961406	!ece/hb
      omfmax=om
      nf=10000
      df=(omfmax-omfmin)/nf
      omf(1)=1e-12
      DO i=2,nf
         omf(i)=(i-1)*df
      ENDDO
      CALL four(cohccint_x,tme,dt,nt,cohccf_x,omf,df,nf)
      CALL four(cohccint_y,tme,dt,nt,cohccf_y,omf,df,nf)
      Open(unit=3,file='cohccft.txt')
      DO iw=1,nf
         WRITE(3,*)omf(iw),cabs(cohccf_x(iw)),cabs(cohccf_y(iw))
      ENDDO
      CLOSE(3)

      CALL four(cohvvint_x,tme,dt,nt,cohvvf_x,omf,df,nf)
      CALL four(cohvvint_y,tme,dt,nt,cohvvf_y,omf,df,nf)
      Open(unit=3,file='cohvvft.txt')
      DO iw=1,nf
         WRITE(3,*)omf(iw),cabs(cohvvf_x(iw)),cabs(cohvvf_y(iw))
      ENDDO
      CLOSE(3)

      CALL four(cohcvint_x,tme,dt,nt,cohcvf_x,omf,df,nf)
      CALL four(cohcvint_y,tme,dt,nt,cohcvf_y,omf,df,nf)
      Open(unit=3,file='cohcvft.txt')
      DO iw=1,nf
         WRITE(3,*)omf(iw),cabs(cohcvf_x(iw)),cabs(cohcvf_y(iw))
      ENDDO
      CLOSE(3)

      ntf=5000
      dtf=(tmax-tmin)/ntf
      DO it=1,ntf
         tmef(it)=tmin+(it-1)*dtf
      ENDDO
      CALL ifour(cohccf_x,omf,df,nf,cohccif_x,tmef,dtf,ntf)
      CALL ifour(cohccf_y,omf,df,nf,cohccif_y,tmef,dtf,ntf)
      Open(unit=4,file='cohccift.txt')
      Do it=1,ntf
      Write(4,*)tmef(it),cohccif_x(it),cohccif_y(it)
      enddo
      Close(4)

      CALL ifour(cohvvf_x,omf,df,nf,cohvvif_x,tmef,dtf,ntf)
      CALL ifour(cohvvf_y,omf,df,nf,cohvvif_y,tmef,dtf,ntf)
      Open(unit=4,file='cohvvift.txt')
      Do it=1,ntf
      Write(4,*)tmef(it),cohvvif_x(it),cohvvif_y(it)
      enddo
      Close(4)

      CALL ifour(cohcvf_x,omf,df,nf,cohcvif_x,tmef,dtf,ntf)
      CALL ifour(cohcvf_y,omf,df,nf,cohcvif_y,tmef,dtf,ntf)
      Open(unit=4,file='cohcvift.txt')
      Do it=1,ntf
      Write(4,*)tmef(it),cohcvif_x(it),cohcvif_y(it)
      enddo
      Close(4)
      END!end photocurre
! bien doi fourier tu time->frequency

      SUBROUTINE four(datat,tme,dt,nt,dataf,omf,df,nf)
      REAL datat(nt), tme(nt), omf(nt)
      COMPLEX dataf(nt), z ,sum1
      z=cmplx(0.0,1.0)
      DO i=1,nf
         sum1=(0.,0.)
         DO j=1,nt
            sum1=sum1+datat(j)*exp(z*omf(i)*tme(j))
         ENDDO
         dataf(i)=sum1*dt
      ENDDO
      END
!bien doi nguoc lai
      SUBROUTINE ifour(dataf,omf,df,nf,datat,tme,dt,nt)
      PARAMETER (pi=3.14159)
      REAL datat(nt),tme(nt),omf(nf)
      COMPLEX dataf(nf),sum1,z
      z=cmplx(0.,1.)
      DO i=1,nt
         sum1=(0.,0.)
         DO j=1,nf
            sum1=sum1+dataf(j)*exp(-z*omf(j)*tme(i))
         ENDDO
         datat(i)=2*real(sum1)*df/(2.*pi)
      ENDDO
      END
! density electron and hole
      SUBROUTINE density(ye,yh,eden,hden)
      PARAMETER (nkx=35,pi=3.14159,NC=4,NV=8)
      REAL sk(nkx),ye(NC,nkx,nkx),yh(NV,nkx,nkx)
      COMMON/block07/nk,sk,dk
      sum1=0.
      DO L=1,NC
         DO i=1,nk
         DO j=1,nk
            sum1=sum1+(dk**2)*ye(L,i,j)
         ENDDO
         ENDDO
      ENDDO
      eden=sum1/((2*pi)**2)
      sum2=0.
      DO L=1,NV
         DO i=1,nk
         DO j=1,nk
            sum2=sum2+(dk**2)*yh(L,i,j)
         ENDDO
         ENDDO
      ENDDO
      hden=sum2/((2*pi)**2)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SBE(t,om,ye,yh,yc,yv,yp,dedt,dhdt,dcdt,dvdt,dpdt)
      PARAMETER(nkx=35,NC=4,NV=8)
      REAL sk(nkx),bandc(NC,nkx,nkx),bandv(NV,nkx,nkx),
     &     ye(NC,nkx,nkx),yh(NV,nkx,nkx),
     &     dedt(NC,nkx,nkx),dhdt(NV,nkx,nkx),
     &     yefd(NC,nkx,nkx),yhfd(NV,nkx,nkx)
      COMPLEX z,occ(NC,NC,nkx,nkx),ovv(NV,NV,nkx,nkx),
     &        ocv(NC,NV,nkx,nkx),ovc(NV,NC,nkx,nkx),
     &        yc(NC,NC,nkx,nkx),yv(NV,NV,nkx,nkx),yp(NC,NV,nkx,nkx),
     &        dcdt(NC,NC,nkx,nkx),dvdt(NV,NV,nkx,nkx),
     &        dpdt(NC,NV,nkx,nkx),sum1
      COMMON/block01/hb,xm0,beta
      COMMON/block07/nk,sk,dk
      COMMON/block08/bandc,bandv
      COMMON/block11/occ,ovv,ocv,ovc

      z=cmplx(0.,1.)
      tau1=200.   !fs
      tau2=200.   !fs
      CALL rabifreq(t,om)
      DO i=1,nk
      DO j=1,nk
         DO ic=1,NC
         DO jc=ic,NC
            sum1=0.
            DO kc=1,NC
               sum1=sum1+occ(kc,ic,i,j)*yc(kc,jc,i,j)
     &                  -occ(jc,kc,i,j)*yc(ic,kc,i,j)
            ENDDO
            DO kv=1,NV
               sum1=sum1+ovc(kv,ic,i,j)*conjg(yp(jc,kv,i,j))
     &                  -ocv(jc,kv,i,j)*yp(ic,kv,i,j)
            ENDDO
            IF (ic.eq.jc) THEN
               dedt(ic,i,j)=z/hb*sum1
               dcdt(ic,ic,i,j)=dedt(ic,i,j)
            ELSE
               dcdt(ic,jc,i,j)=z/hb*sum1
     &              +z/hb*(bandc(ic,i,j)-bandc(jc,i,j))*yc(ic,jc,i,j)
     &              -yc(ic,jc,i,j)/tau2
               dcdt(jc,ic,i,j)=conjg(dcdt(ic,jc,i,j))
            ENDIF
         ENDDO
         ENDDO
         DO iv=1,NV
         DO jv=iv,NV
            sum1=0.
            DO kc=1,NC
               sum1=sum1+ocv(kc,iv,i,j)*yp(kc,jv,i,j)
     &                  -ovc(jv,kc,i,j)*conjg(yp(kc,iv,i,j))
            ENDDO
            DO kv=1,NV
               sum1=sum1+ovv(kv,iv,i,j)*yv(kv,jv,i,j)
     &                  -ovv(jv,kv,i,j)*yv(iv,kv,i,j)
            ENDDO
            IF (iv.eq.jv) THEN
               dhdt(iv,i,j)=-z/hb*sum1
               dvdt(iv,iv,i,j)=-dhdt(iv,i,j)
            ELSE
               dvdt(iv,jv,i,j)=z/hb*sum1
     &              +z/hb*(bandv(iv,i,j)-bandv(jv,i,j))*yv(iv,jv,i,j)
     &              -yv(iv,jv,i,j)/tau2
               dvdt(jv,iv,i,j)=conjg(dvdt(iv,jv,i,j))
            ENDIF
         ENDDO
         ENDDO
         DO ic=1,NC
         DO iv=1,NV
            sum1=0.
            DO kc=1,NC
               sum1=sum1+occ(kc,ic,i,j)*yp(kc,iv,i,j)
     &                  -ovc(iv,kc,i,j)*yc(ic,kc,i,j)
            ENDDO
            DO kv=1,NV
               sum1=sum1+ovc(kv,ic,i,j)*yv(kv,iv,i,j)
     &                  -ovv(iv,kv,i,j)*yp(ic,kv,i,j)
            ENDDO
            dpdt(ic,iv,i,j)=z/hb*sum1
     &           +z/hb*(bandc(ic,i,j)-bandv(iv,i,j))*yp(ic,iv,i,j)
     &           -yp(ic,iv,i,j)/tau2
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      END! end SBE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RK2SBE(t,h,om,ye,yh,yc,yv,yp,
     &                  yeout,yhout,ycout,yvout,ypout)
      PARAMETER(nkx=35,NC=4,NV=8)
      REAL sk(nkx),ye(NC,nkx,nkx),yeout(NC,nkx,nkx),dedt(NC,nkx,nkx),
     &     yet(NC,nkx,nkx),yh(NV,nkx,nkx),yhout(NV,nkx,nkx),
     &     dhdt(NV,nkx,nkx),yht(NV,nkx,nkx)
      COMPLEX yc(NC,NC,nkx,nkx),dcdt(NC,NC,nkx,nkx),
     &        ycout(NC,NC,nkx,nkx),yct(NC,NC,nkx,nkx),
     &        yv(NV,NV,nkx,nkx),dvdt(NV,NV,nkx,nkx),
     &        yvout(NV,NV,nkx,nkx),yvt(NV,NV,nkx,nkx),
     &        yp(NC,NV,nkx,nkx),dpdt(NC,NV,nkx,nkx),
     &        ypout(NC,NV,nkx,nkx),ypt(NC,NV,nkx,nkx)
      COMMON/block07/nk,sk,dk
      CALL SBE(t,om,ye,yh,yc,yv,yp,dedt,dhdt,dcdt,dvdt,dpdt)
      hh=h*0.5
      th=t+hh
      DO i=1,nk
      DO j=1,nk
         DO L=1,NC
            yet(L,i,j)=ye(L,i,j)+hh*dedt(L,i,j)
            yct(L,L,i,j)=yet(L,i,j)
         ENDDO
         DO L1=1,NC-1
         DO L2=L1+1,NC
            yct(L1,L2,i,j)=yc(L1,L2,i,j)+hh*dcdt(L1,L2,i,j)
            yct(L2,L1,i,j)=conjg(yct(L1,L2,i,j))
         ENDDO
         ENDDO
         DO L=1,NV
            yht(L,i,j)=yh(L,i,j)+hh*dhdt(L,i,j)
            yvt(L,L,i,j)=1.-yht(L,i,j)
         ENDDO
         DO L1=1,NV-1
         DO L2=L1+1,NV
            yvt(L1,L2,i,j)=yv(L1,L2,i,j)+hh*dvdt(L1,L2,i,j)
            yvt(L2,L1,i,j)=conjg(yvt(L1,L2,i,j))
         ENDDO
         ENDDO
         DO Lc=1,NC
         DO Lv=1,NV
            ypt(Lc,Lv,i,j)=yp(Lc,Lv,i,j)+hh*dpdt(Lc,Lv,i,j)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      CALL SBE(th,om,yet,yht,yct,yvt,ypt,dedt,dhdt,dcdt,dvdt,dpdt)
      DO i=1,nk
      DO j=1,nk
         DO L=1,NC
            yeout(L,i,j)=ye(L,i,j)+h*dedt(L,i,j)
            ycout(L,L,i,j)=yeout(L,i,j)
         ENDDO
         DO L1=1,NC-1
         DO L2=L1+1,NC
            ycout(L1,L2,i,j)=yc(L1,L2,i,j)+h*dcdt(L1,L2,i,j)
            ycout(L2,L1,i,j)=conjg(ycout(L1,L2,i,j))
         ENDDO
         ENDDO
         DO L=1,NV
            yhout(L,i,j)=yh(L,i,j)+h*dhdt(L,i,j)
            yvout(L,L,i,j)=1.-yhout(L,i,j)
         ENDDO
         DO L1=1,NV-1
         DO L2=L1+1,NV
            yvout(L1,L2,i,j)=yv(L1,L2,i,j)+h*dvdt(L1,L2,i,j)
            yvout(L2,L1,i,j)=conjg(yvout(L1,L2,i,j))
         ENDDO
         ENDDO
         DO Lc=1,NC
         DO Lv=1,NV
            ypout(Lc,Lv,i,j)=yp(Lc,Lv,i,j)+h*dpdt(Lc,Lv,i,j)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      END! end rk2sbe
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE bandstr(HxxxBE,Hxxx,Pxxx)
      PARAMETER(NBX=22,nkx=35,NC=4,NV=8)
      REAL sk(nkx),val(14*NBX),WORK(14*14*NBX),
     &     pot8c(NBX,NBX),pot7c(NBX,NBX),pot6c(NBX,NBX),
     &     pot8v(NBX,NBX),pot7v(NBX,NBX),dmm(NBX,NBX),p0m(NBX,NBX),
     &     p1m(NBX,NBX),qm(NBX,NBX),gcm(NBX,NBX),g1m(NBX,NBX),
     &     g2m(NBX,NBX),g3m(NBX,NBX),ckm(NBX,NBX),
     &     bandc(NC,nkx,nkx),bandv(NV,nkx,nkx)
      COMPLEX kz(NBX,NBX),ham0(14*NBX,14*NBX),ham(14*NBX,14*NBX),
     &        vec(14*NBX,14*NBX),Px(14*NBX,14*NBX),Py(14*NBX,14*NBX),
     &        Pz(14*NBX,14*NBX),sum1,sum2,sum3,pcv(NC,NV,3,nkx,nkx),
     &        pcc(NC,NC,3,nkx,nkx),pvv(NV,NV,3,nkx,nkx)
      EXTERNAL HxxxBE,Hxxx,Pxxx
      COMMON/block01/hb,xm0,beta
      COMMON/block07/nk,sk,dk
      COMMON/block08/bandc,bandv
      COMMON/block09/pcc,pvv,pcv

      NB=NBX    !the number of plane waves
      CALL BME(NB,kz,pot8c,pot7c,pot6c,pot8v,pot7v,
     &         dmm,p0m,p1m,qm,gcm,g1m,g2m,g3m,ckm)
      CALL HxxxBE(NB,kz,pot8c,pot7c,pot6c,pot8v,pot7v,
     &            dmm,p0m,p1m,qm,gcm,g1m,g2m,g3m,ckm,ham0)
      DO i=1,nk
         akx=sk(i)
         DO j=1,nk
            aky=sk(j)
            CALL Hxxx(NB,akx,aky,kz,dmm,p0m,p1m,qm,gcm,g1m,g2m,g3m,
     &                ckm,ham0,ham)
            CALL CHIEV(ham,14*NB,14*NB,val,vec,14*NB,WORK,1,INFO)
            DO Lc=1,NC
               bandc(Lc,i,j)=val(6*NB+Lc)
            ENDDO
            DO Lv=1,NV
               bandv(Lv,i,j)=val(6*NB+1-Lv)
            ENDDO
            CALL Pxxx(NB,akx,aky,kz,dmm,p0m,p1m,qm,gcm,g1m,g2m,g3m,ckm,
     &                Px,Py,Pz)
            DO Lc=1,NC
               DO Lv=1,NV
                  sum1=(0.,0.)
                  sum2=(0.,0.)
                  sum3=(0.,0.)
                  DO n=1,14*NB
                  DO m=1,14*NB
                     sum1=sum1
     &               +conjg(vec(m,6*NB+Lc))*Px(m,n)*vec(n,6*NB+1-Lv)
                     sum2=sum2
     &               +conjg(vec(m,6*NB+Lc))*Py(m,n)*vec(n,6*NB+1-Lv)
                     sum3=sum3
     &               +conjg(vec(m,6*NB+Lc))*Pz(m,n)*vec(n,6*NB+1-Lv)
                  ENDDO
                  ENDDO
                  pcv(Lc,Lv,1,i,j)=sum1*xm0/hb
                  pcv(Lc,Lv,2,i,j)=sum2*xm0/hb
                  pcv(Lc,Lv,3,i,j)=sum3*xm0/hb
               ENDDO
            ENDDO
            DO L1=1,NC
               DO L2=L1,NC
                  sum1=(0.,0.)
                  sum2=(0.,0.)
                  sum3=(0.,0.)
                  DO n=1,14*NB
                  DO m=1,14*NB
                     sum1=sum1
     &               +conjg(vec(m,6*NB+L1))*Px(m,n)*vec(n,6*NB+L2)
                     sum2=sum2
     &               +conjg(vec(m,6*NB+L1))*Py(m,n)*vec(n,6*NB+L2)
                     sum3=sum3
     &               +conjg(vec(m,6*NB+L1))*Pz(m,n)*vec(n,6*NB+L2)
                  ENDDO
                  ENDDO
                  pcc(L1,L2,1,i,j)=sum1*xm0/hb
                  pcc(L2,L1,1,i,j)=conjg(sum1)*xm0/hb
                  pcc(L1,L2,2,i,j)=sum2*xm0/hb
                  pcc(L2,L1,2,i,j)=conjg(sum2)*xm0/hb
                  pcc(L1,L2,3,i,j)=sum3*xm0/hb
                  pcc(L2,L1,3,i,j)=conjg(sum3)*xm0/hb
               ENDDO
            ENDDO
            DO L1=1,NV
               DO L2=L1,NV
                  sum1=(0.,0.)
                  sum2=(0.,0.)
                  sum3=(0.,0.)
                  DO n=1,14*NB
                  DO m=1,14*NB
                     sum1=sum1
     &               +conjg(vec(m,6*NB+1-L1))*Px(m,n)*vec(n,6*NB+1-L2)
                     sum2=sum2
     &               +conjg(vec(m,6*NB+1-L1))*Py(m,n)*vec(n,6*NB+1-L2)
                     sum3=sum3
     &               +conjg(vec(m,6*NB+1-L1))*Pz(m,n)*vec(n,6*NB+1-L2)
                  ENDDO
                  ENDDO
                  pvv(L1,L2,1,i,j)=sum1*xm0/hb
                  pvv(L2,L1,1,i,j)=conjg(sum1)*xm0/hb
                  pvv(L1,L2,2,i,j)=sum2*xm0/hb
                  pvv(L2,L1,2,i,j)=conjg(sum2)*xm0/hb
                  pvv(L1,L2,3,i,j)=sum3*xm0/hb
                  pvv(L2,L1,3,i,j)=conjg(sum3)*xm0/hb
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      OPEN(UNIT=1,FILE='band110qw.txt')
      DO Lc=1,NC
         DO i=1,(nk+1)/2
            WRITE(1,*)sk(i),bandc(Lc,i,(nk+1)/2)
         ENDDO
         DO j=(nk+1)/2+1,nk
            WRITE(1,*)sk(j),bandc(Lc,(nk+1)/2,j)
         ENDDO
         WRITE(1,*)' '
      ENDDO
      DO Lv=1,NV
         DO i=1,(nk+1)/2
            WRITE(1,*)sk(i),bandv(Lv,i,(nk+1)/2)
         ENDDO
         DO j=(nk+1)/2+1,nk
            WRITE(1,*)sk(j),bandv(Lv,(nk+1)/2,j)
         ENDDO
         WRITE(1,*)' '
      ENDDO
      CLOSE(1)

      OPEN(UNIT=2,FILE='momentum.txt')
      DO i=1,nk
         DO j=1,nk
            WRITE(2,'(22f23.16)')sk(i),sk(j),
     &              cabs(pcv(1,1,1,i,j))+cabs(pcv(1,2,1,i,j))
     &             +cabs(pcv(2,1,1,i,j))+cabs(pcv(2,2,1,i,j))
     &             ,cabs(pcv(1,3,1,i,j))+cabs(pcv(1,4,1,i,j))
     &             +cabs(pcv(2,3,1,i,j))+cabs(pcv(2,4,1,i,j))
     &             ,cabs(pcv(1,1,2,i,j))+cabs(pcv(1,2,2,i,j))
     &             +cabs(pcv(2,1,2,i,j))+cabs(pcv(2,2,2,i,j))
     &             ,cabs(pcv(1,3,2,i,j))+cabs(pcv(1,4,2,i,j))
     &             +cabs(pcv(2,3,2,i,j))+cabs(pcv(2,4,2,i,j))
         ENDDO
         WRITE(2,*)' '
      ENDDO
      CLOSE(2)
      END!end bandstruc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE H110BE(NB,kz,pot8c,pot7c,pot6c,pot8v,pot7v,
     &                  dmm,p0m,p1m,qm,gcm,g1m,g2m,g3m,ckm,ham)
      REAL pot8c(NB,NB),pot7c(NB,NB),pot6c(NB,NB),pot8v(NB,NB),
     &     pot7v(NB,NB),dmm(NB,NB),p0m(NB,NB),p1m(NB,NB),qm(NB,NB),
     &     gcm(NB,NB),g1m(NB,NB),g2m(NB,NB),g3m(NB,NB),ckm(NB,NB)
      COMPLEX z,kz(NB,NB),ham(14*NB,14*NB),sum4,sum5,sum6,sum7
      COMMON/block01/hb,xm0,beta
      z=cmplx(0.,1.)
      DO i=1,14*NB
         DO j=1,14*NB
            ham(i,j)=(0.,0.)
         ENDDO
      ENDDO
      DO i=1,NB
         DO j=1,NB
! H8c8c
            ham(i,j)=pot8c(i,j)
            ham(i+NB,j+NB)=pot8c(i,j)
            ham(i+2*NB,j+2*NB)=pot8c(i,j)
            ham(i+3*NB,j+3*NB)=pot8c(i,j)
! H7c7c
            ham(i+4*NB,j+4*NB)=pot7c(i,j)
            ham(i+5*NB,j+5*NB)=pot7c(i,j)
! H6c6c
            sum1=0.
            DO k=1,NB
               DO l=1,NB
                  sum1=sum1+kz(i,k)*gcm(k,l)*kz(l,j)
               ENDDO
            ENDDO
            ham(i+6*NB,j+6*NB)=beta*sum1+pot6c(i,j)
            ham(i+7*NB,j+7*NB)=beta*sum1+pot6c(i,j)
! H8v8v
            sum1=0.
            sum2=0.
            sum3=0.
            DO k=1,NB
               DO l=1,NB
                  sum1=sum1+kz(i,k)
     &                *((2*g1m(k,l)-g2m(k,l)-3*g3m(k,l))/2)*kz(l,j)
                  sum2=sum2+kz(i,k)
     &                *((2*g1m(k,l)+g2m(k,l)+3*g3m(k,l))/2)*kz(l,j)
                  sum3=sum3+kz(i,k)*(g2m(k,l)-g3m(k,l))*kz(l,j)
               ENDDO
            ENDDO
            sum4=(0.,0.)
            DO l=1,NB
               sum4=sum4+.5*(kz(i,l)*ckm(l,j)+ckm(i,l)*kz(l,j))
            ENDDO
            ham(i+8*NB,j+8*NB)=-beta*sum1+pot8v(i,j)
            ham(i+9*NB,j+9*NB)=-beta*sum2+pot8v(i,j)
            ham(i+10*NB,j+10*NB)=-beta*sum2+pot8v(i,j)
            ham(i+11*NB,j+11*NB)=-beta*sum1+pot8v(i,j)
            ham(i+8*NB,j+10*NB)=-.5*sqrt(3.)*beta*sum3
            ham(i+9*NB,j+11*NB)=-.5*sqrt(3.)*beta*sum3
            ham(i+8*NB,j+9*NB)=-.25*z*sum4
            ham(i+8*NB,j+11*NB)=-3*sqrt(3.)/4*z*sum4
            ham(i+9*NB,j+10*NB)=sqrt(3.)/4*z*sum4
            ham(i+10*NB,j+11*NB)=-.25*z*sum4
! H7v7v
            sum1=0.
            DO k=1,NB
               DO l=1,NB
                  sum1=sum1+kz(i,k)*g1m(k,l)*kz(l,j)
               ENDDO
            ENDDO
            ham(i+12*NB,j+12*NB)=-beta*sum1+pot7v(i,j)
            ham(i+13*NB,j+13*NB)=-beta*sum1+pot7v(i,j)
! H8c6c
            sum5=(0.,0.)
            DO l=1,NB
               sum5=sum5+.5*(kz(i,l)*p1m(l,j)+p1m(i,l)*kz(l,j))
            ENDDO
            ham(i+NB,j+6*NB)=-sqrt(2./3)*z*sum5
            ham(i+2*NB,j+7*NB)=-sqrt(2./3)*z*sum5
! H8c8v
            sum6=(0.,0.)
            DO l=1,NB
               sum6=sum6+.5*(kz(i,l)*qm(l,j)+qm(i,l)*kz(l,j))
            ENDDO
            ham(i,j+9*NB)=sqrt(1./3)*sum6
            ham(i+NB,j+8*NB)=sqrt(1./3)*sum6
            ham(i+2*NB,j+11*NB)=-sqrt(1./3)*sum6
            ham(i+3*NB,j+10*NB)=-sqrt(1./3)*sum6
! H8c7v
            ham(i,j+12*NB)=-sqrt(1./6)*sum6
            ham(i+NB,j+13*NB)=sqrt(.5)*sum6
            ham(i+2*NB,j+12*NB)=sqrt(.5)*sum6
            ham(i+3*NB,j+13*NB)=-sqrt(1./6)*sum6
! H7c6c
            ham(i+4*NB,j+6*NB)=sqrt(1./3)*z*sum5
            ham(i+5*NB,j+7*NB)=-sqrt(1./3)*z*sum5
! H7c8v
            ham(i+4*NB,j+8*NB)=-sqrt(1./6)*sum6
            ham(i+4*NB,j+10*NB)=sqrt(.5)*sum6
            ham(i+5*NB,j+9*NB)=sqrt(.5)*sum6
            ham(i+5*NB,j+11*NB)=-sqrt(1./6)*sum6
! H6c8v
            sum7=(0.,0.)
            DO l=1,NB
               sum7=sum7+.5*(kz(i,l)*p0m(l,j)+p0m(i,l)*kz(l,j))
            ENDDO
            ham(i+6*NB,j+9*NB)=sqrt(2./3)*sum7
            ham(i+7*NB,j+10*NB)=sqrt(2./3)*sum7
! H6c7v
            ham(i+6*NB,j+12*NB)=-sqrt(1./3)*sum7
            ham(i+7*NB,j+13*NB)=sqrt(1./3)*sum7
! H8v7v
            sum1=0.
            sum2=0.
            DO k=1,NB
               DO l=1,NB
                  sum1=sum1+kz(i,k)*(g2m(k,l)-g3m(k,l))*kz(l,j)
                  sum2=sum2+kz(i,k)*(g2m(k,l)+3*g3m(k,l))*kz(l,j)
               ENDDO
            ENDDO
            ham(i+8*NB,j+12*NB)=-.5*sqrt(.5)*z*sum4
            ham(i+8*NB,j+13*NB)=sqrt(3./2)*beta*sum1
            ham(i+9*NB,j+12*NB)=sqrt(.5)*beta*sum2
            ham(i+9*NB,j+13*NB)=.5*sqrt(3./2)*z*sum4
            ham(i+10*NB,j+12*NB)=.5*sqrt(3./2)*z*sum4
            ham(i+10*NB,j+13*NB)=-sqrt(.5)*beta*sum2
            ham(i+11*NB,j+12*NB)=-sqrt(3./2)*beta*sum1
            ham(i+11*NB,j+13*NB)=-.5*sqrt(.5)*z*sum4
         ENDDO
      ENDDO
      DO i=1,14*NB-1
         DO j=i+1,14*NB
            ham(j,i)=conjg(ham(i,j))
         ENDDO
      ENDDO
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE H110(NB,kx,ky,kz,dmm,p0m,p1m,qm,gcm,g1m,g2m,g3m,
     &                 ckm,ham0,ham)
      REAL kx,ky,dmm(NB,NB),p0m(NB,NB),p1m(NB,NB),qm(NB,NB),
     &     gcm(NB,NB),g1m(NB,NB),g2m(NB,NB),g3m(NB,NB),ckm(NB,NB)
      COMPLEX z,kz(NB,NB),kp,km,sum2,sum3,
     &        ham0(14*NB,14*NB),ham(14*NB,14*NB)
      COMMON/block01/hb,xm0,beta
      z=cmplx(0.,1.)
      kp=kx+z*ky
      km=kx-z*ky
      DO i=1,14*NB
         DO j=1,14*NB
            ham(i,j)=(0.,0)
         ENDDO
      ENDDO
      DO i=1,NB
         DO j=1,NB
! H8c6c
            ham(i,j+6*NB)=sqrt(.5)*z*p1m(i,j)*km
            ham(i+NB,j+7*NB)=sqrt(1./6)*z*p1m(i,j)*km
            ham(i+2*NB,j+6*NB)=-sqrt(1./6)*z*p1m(i,j)*kp
            ham(i+3*NB,j+7*NB)=-sqrt(.5)*z*p1m(i,j)*kp
! H8c8v
            ham(i,j+8*NB)=z*dmm(i,j)/3+.5*qm(i,j)*kx
            ham(i,j+10*NB)=.5/sqrt(3.)*qm(i,j)*(kx+2*z*ky)
            ham(i+NB,j+9*NB)=z*dmm(i,j)/3-.5*qm(i,j)*kx
            ham(i+NB,j+11*NB)=.5/sqrt(3.)*qm(i,j)*(kx+2*z*ky)
            ham(i+2*NB,j+8*NB)=.5/sqrt(3.)*qm(i,j)*(kx-2*z*ky)
            ham(i+2*NB,j+10*NB)=z*dmm(i,j)/3-.5*qm(i,j)*kx
            ham(i+3*NB,j+9*NB)=.5/sqrt(3.)*qm(i,j)*(kx-2*z*ky)
            ham(i+3*NB,j+11*NB)=z*dmm(i,j)/3+.5*qm(i,j)*kx
! H8c7v
            ham(i,j+13*NB)=-sqrt(1./6)*qm(i,j)*(kx+2*z*ky)
            ham(i+NB,j+12*NB)=sqrt(.5)*qm(i,j)*kx
            ham(i+2*NB,j+13*NB)=-sqrt(.5)*qm(i,j)*kx
            ham(i+3*NB,j+12*NB)=sqrt(1./6)*qm(i,j)*(kx-2*z*ky)
! H7c6c
            ham(i+4*NB,j+7*NB)=sqrt(1./3)*z*p1m(i,j)*km
            ham(i+5*NB,j+6*NB)=sqrt(1./3)*z*p1m(i,j)*kp
! H7c8v
            ham(i+4*NB,j+9*NB)=sqrt(.5)*qm(i,j)*kx
            ham(i+4*NB,j+11*NB)=sqrt(1./6)*qm(i,j)*(kx+2*z*ky)
            ham(i+5*NB,j+8*NB)=-sqrt(1./6)*qm(i,j)*(kx-2*z*ky)
            ham(i+5*NB,j+10*NB)=-sqrt(.5)*qm(i,j)*kx
! H7c7v
            ham(i+4*NB,j+12*NB)=-2*z*dmm(i,j)/3
            ham(i+5*NB,j+13*NB)=-2*z*dmm(i,j)/3
! H6c6c
            ham(i+6*NB,j+6*NB)=beta*gcm(i,j)*(kx**2+ky**2)
            ham(i+7*NB,j+7*NB)=beta*gcm(i,j)*(kx**2+ky**2)
! H6c8v
            ham(i+6*NB,j+8*NB)=-sqrt(.5)*p0m(i,j)*kp
            ham(i+6*NB,j+10*NB)=sqrt(1./6)*p0m(i,j)*km
            ham(i+7*NB,j+9*NB)=-sqrt(1./6)*p0m(i,j)*kp
            ham(i+7*NB,j+11*NB)=sqrt(.5)*p0m(i,j)*km
! H6c7v
            ham(i+6*NB,j+13*NB)=-sqrt(1./3)*p0m(i,j)*km
            ham(i+7*NB,j+12*NB)=-sqrt(1./3)*p0m(i,j)*kp
! H8v8v
            sum2=(0.,0.)
            sum3=(0.,0.)
            DO l=1,NB
               sum2=sum2+.5*(kz(i,l)*g2m(l,j)+g2m(i,l)*kz(l,j))
               sum3=sum3+.5*(kz(i,l)*g3m(l,j)+g3m(i,l)*kz(l,j))
            ENDDO
            ham(i+8*NB,j+8*NB)=-beta*(g1m(i,j)+g2m(i,j))*kx**2
     &              -beta*((2*g1m(i,j)-g2m(i,j)+3*g3m(i,j))/2)*ky**2
     &                         -sqrt(3.)/4*ckm(i,j)*ky
            ham(i+8*NB,j+9*NB)=2*sqrt(3.)*beta*(sum3*kx-z*sum2*ky)
            ham(i+8*NB,j+10*NB)=sqrt(3.)*beta*g2m(i,j)*kx**2
     &                   -.5*sqrt(3.)*beta*(g2m(i,j)+g3m(i,j))*ky**2
     &                   -2*sqrt(3.)*z*beta*g3m(i,j)*kx*ky
     &                   +.25*ckm(i,j)*(ky+4*z*kx)
            ham(i+9*NB,j+9*NB)=-beta*(g1m(i,j)-g2m(i,j))*kx**2
     &              -beta*((2*g1m(i,j)+g2m(i,j)-3*g3m(i,j))/2)*ky**2
     &                         +3*sqrt(3.)/4*ckm(i,j)*ky
            ham(i+9*NB,j+11*NB)=sqrt(3.)*beta*g2m(i,j)*kx**2
     &                   -.5*sqrt(3.)*beta*(g2m(i,j)+g3m(i,j))*ky**2
     &                   -2*sqrt(3.)*z*beta*g3m(i,j)*kx*ky
     &                   -.25*ckm(i,j)*(ky+4*z*kx)
            ham(i+10*NB,j+10*NB)=-beta*(g1m(i,j)-g2m(i,j))*kx**2
     &              -beta*((2*g1m(i,j)+g2m(i,j)-3*g3m(i,j))/2)*ky**2
     &                           -3*sqrt(3.)/4*ckm(i,j)*ky
            ham(i+10*NB,j+11*NB)=-2*sqrt(3.)*beta*(sum3*kx-z*sum2*ky)
            ham(i+11*NB,j+11*NB)=-beta*(g1m(i,j)+g2m(i,j))*kx**2
     &              -beta*((2*g1m(i,j)-g2m(i,j)+3*g3m(i,j))/2)*ky**2
     &                           +sqrt(3.)/4*ckm(i,j)*ky
! H8v7v
            ham(i+8*NB,j+12*NB)=-sqrt(6.)*beta*(sum3*kx-z*sum2*ky)
            ham(i+8*NB,j+13*NB)=-sqrt(6.)*beta*g2m(i,j)*kx**2
     &                    +sqrt(3./2)*beta*(g2m(i,j)+g3m(i,j))*ky**2
     &                    +2*sqrt(6.)*z*beta*g3m(i,j)*kx*ky
     &                    +.5*sqrt(.5)*ckm(i,j)*(2*ky-z*kx)
            ham(i+9*NB,j+12*NB)=-sqrt(2.)*beta*g2m(i,j)*kx**2
     &                    +sqrt(.5)*beta*(g2m(i,j)-3*g3m(i,j))*ky**2
     &                    +.5*sqrt(3./2)*z*ckm(i,j)*kx
            ham(i+9*NB,j+13*NB)=3*sqrt(2.)*beta*(sum3*kx-z*sum2*ky)
            ham(i+10*NB,j+12*NB)=3*sqrt(2.)*beta*(sum3*kx+z*sum2*ky)
            ham(i+10*NB,j+13*NB)=sqrt(2.)*beta*g2m(i,j)*kx**2
     &                    -sqrt(.5)*beta*(g2m(i,j)-3*g3m(i,j))*ky**2
     &                    -.5*sqrt(3./2)*z*ckm(i,j)*kx
            ham(i+11*NB,j+12*NB)=sqrt(6.)*beta*g2m(i,j)*kx**2
     &                    -sqrt(3./2)*beta*(g2m(i,j)+g3m(i,j))*ky**2
     &                    +2*sqrt(6.)*z*beta*g3m(i,j)*kx*ky
     &                    +.5*sqrt(.5)*ckm(i,j)*(2*ky+z*kx)
            ham(i+11*NB,j+13*NB)=-sqrt(6.)*beta*(sum3*kx+z*sum2*ky)
! H7v7v
            ham(i+12*NB,j+12*NB)=-beta*g1m(i,j)*(kx**2+ky**2)
            ham(i+13*NB,j+13*NB)=-beta*g1m(i,j)*(kx**2+ky**2)
         ENDDO
      ENDDO
      DO i=1,14*NB-1
         DO j=i+1,14*NB
            ham(j,i)=conjg(ham(i,j))
         ENDDO
      ENDDO
      DO i=1,14*NB
         DO j=1,14*NB
            ham(i,j)=ham0(i,j)+ham(i,j)
         ENDDO
      ENDDO
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE P110(NB,kx,ky,kz,dmm,p0m,p1m,qm,gcm,g1m,g2m,g3m,ckm,
     &                Px,Py,Pz)
      REAL kx,ky,dmm(NB,NB),p0m(NB,NB),p1m(NB,NB),qm(NB,NB),
     &     gcm(NB,NB),g1m(NB,NB),g2m(NB,NB),g3m(NB,NB),
     &     ckm(NB,NB)
      COMPLEX z,kz(NB,NB),sum1,sum2,sum3,sum4,
     &        Px(14*NB,14*NB),Py(14*NB,14*NB),Pz(14*NB,14*NB)
      COMMON/block01/hb,xm0,beta
      z=cmplx(0.,1.)
      DO i=1,14*NB
         DO j=1,14*NB
            Px(i,j)=(0.,0.)
            Py(i,j)=(0.,0.)
            Pz(i,j)=(0.,0.)
         ENDDO
      ENDDO
      DO i=1,NB
         DO j=1,NB
! 8c6c
            Px(i,j+6*NB)=sqrt(1./2)*z*p1m(i,j)
            Px(i+NB,j+7*NB)=sqrt(1./6)*z*p1m(i,j)
            Px(i+2*NB,j+6*NB)=-sqrt(1./6)*z*p1m(i,j)
            Px(i+3*NB,j+7*NB)=-sqrt(1./2)*z*p1m(i,j)
            Py(i,j+6*NB)=sqrt(1./2)*p1m(i,j)
            Py(i+NB,j+7*NB)=sqrt(1./6)*p1m(i,j)
            Py(i+2*NB,j+6*NB)=sqrt(1./6)*p1m(i,j)
            Py(i+3*NB,j+7*NB)=sqrt(1./2)*p1m(i,j)
            Pz(i+NB,j+6*NB)=-sqrt(2./3)*z*p1m(i,j)
            Pz(i+2*NB,j+7*NB)=-sqrt(2./3)*z*p1m(i,j)
! 8c8v
            Px(i,j+8*NB)=.5*qm(i,j)
            Px(i,j+10*NB)=1./(2*sqrt(3.))*qm(i,j)
            Px(i+NB,j+9*NB)=-.5*qm(i,j)
            Px(i+NB,j+11*NB)=1./(2*sqrt(3.))*qm(i,j)
            Px(i+2*NB,j+8*NB)=1./(2*sqrt(3.))*qm(i,j)
            Px(i+2*NB,j+10*NB)=-.5*qm(i,j)
            Px(i+3*NB,j+9*NB)=1./(2*sqrt(3.))*qm(i,j)
            Px(i+3*NB,j+11*NB)=.5*qm(i,j)
            Py(i,j+10*NB)=sqrt(1./3)*z*qm(i,j)
            Py(i+NB,j+11*NB)=sqrt(1./3)*z*qm(i,j)
            Py(i+2*NB,j+8*NB)=-sqrt(1./3)*z*qm(i,j)
            Py(i+3*NB,j+9*NB)=-sqrt(1./3)*z*qm(i,j)
            Pz(i,j+9*NB)=sqrt(1./3)*qm(i,j)
            Pz(i+NB,j+8*NB)=sqrt(1./3)*qm(i,j)
            Pz(i+2*NB,j+11*NB)=-sqrt(1./3)*qm(i,j)
            Pz(i+3*NB,j+10*NB)=-sqrt(1./3)*qm(i,j)
! 8c7v
            Px(i,j+13*NB)=-sqrt(1./6)*qm(i,j)
            Px(i+NB,j+12*NB)=sqrt(1./2)*qm(i,j)
            Px(i+2*NB,j+13*NB)=-sqrt(1./2)*qm(i,j)
            Px(i+3*NB,j+12*NB)=sqrt(1./6)*qm(i,j)
            Py(i,j+13*NB)=-2*sqrt(1./6)*z*qm(i,j)
            Py(i+3*NB,j+12*NB)=-2*sqrt(1./6)*z*qm(i,j)
            Pz(i,j+12*NB)=-sqrt(1./6)*qm(i,j)
            Pz(i+NB,j+13*NB)=sqrt(1./2)*qm(i,j)
            Pz(i+2*NB,j+12*NB)=sqrt(1./2)*qm(i,j)
            Pz(i+3*NB,j+13*NB)=-sqrt(1./6)*qm(i,j)
! 7c6c
            Px(i+4*NB,j+7*NB)=sqrt(1./3)*z*p1m(i,j)
            Px(i+5*NB,j+6*NB)=sqrt(1./3)*z*p1m(i,j)
            Py(i+4*NB,j+7*NB)=sqrt(1./3)*p1m(i,j)
            Py(i+5*NB,j+6*NB)=-sqrt(1./3)*p1m(i,j)
            Pz(i+4*NB,j+6*NB)=sqrt(1./3)*z*p1m(i,j)
            Pz(i+5*NB,j+7*NB)=-sqrt(1./3)*z*p1m(i,j)
! 7c8v
            Px(i+4*NB,j+9*NB)=sqrt(1./2)*qm(i,j)
            Px(i+4*NB,j+11*NB)=sqrt(1./6)*qm(i,j)
            Px(i+5*NB,j+8*NB)=-sqrt(1./6)*qm(i,j)
            Px(i+5*NB,j+10*NB)=-sqrt(1./2)*qm(i,j)
            Py(i+4*NB,j+11*NB)=2*sqrt(1./6)*z*qm(i,j)
            Py(i+5*NB,j+8*NB)=2*sqrt(1./6)*z*qm(i,j)
            Pz(i+4*NB,j+8*NB)=-sqrt(1./6)*qm(i,j)
            Pz(i+4*NB,j+10*NB)=sqrt(1./2)*qm(i,j)
            Pz(i+5*NB,j+9*NB)=sqrt(1./2)*qm(i,j)
            Pz(i+5*NB,j+11*NB)=-sqrt(1./6)*qm(i,j)
! 6c6c
            Px(i+6*NB,i+6*NB)=2*beta*gcm(i,j)*kx
            Px(i+7*NB,i+7*NB)=2*beta*gcm(i,j)*kx
            Py(i+6*NB,i+6*NB)=2*beta*gcm(i,j)*ky
            Py(i+7*NB,i+7*NB)=2*beta*gcm(i,j)*ky
            sum1=(0.,0.)
            DO l=1,NB
               sum1=sum1+.5*(kz(i,l)*gcm(l,j)+gcm(i,l)*kz(l,j))
            ENDDO
            Pz(i+6*NB,j+6*NB)=2*beta*sum1
            Pz(i+7*NB,j+7*NB)=2*beta*sum1
! 6c8v
            Px(i+6*NB,j+8*NB)=-sqrt(1./2)*p0m(i,j)
            Px(i+6*NB,j+10*NB)=sqrt(1./6)*p0m(i,j)
            Px(i+7*NB,j+9*NB)=-sqrt(1./6)*p0m(i,j)
            Px(i+7*NB,j+11*NB)=sqrt(1./2)*p0m(i,j)
            Py(i+6*NB,j+8*NB)=-sqrt(1./2)*z*p0m(i,j)
            Py(i+6*NB,j+10*NB)=-sqrt(1./6)*z*p0m(i,j)
            Py(i+7*NB,j+9*NB)=-sqrt(1./6)*z*p0m(i,j)
            Py(i+7*NB,j+11*NB)=-sqrt(1./2)*z*p0m(i,j)
            Pz(i+6*NB,j+9*NB)=sqrt(2./3)*p0m(i,j)
            Pz(i+7*NB,j+10*NB)=sqrt(2./3)*p0m(i,j)
! 6c7v
            Px(i+6*NB,j+13*NB)=-sqrt(1./3)*p0m(i,j)
            Px(i+7*NB,j+12*NB)=-sqrt(1./3)*p0m(i,j)
            Py(i+6*NB,j+13*NB)=sqrt(1./3)*z*p0m(i,j)
            Py(i+7*NB,j+12*NB)=-sqrt(1./3)*z*p0m(i,j)
            Pz(i+6*NB,j+12*NB)=-sqrt(1./3)*p0m(i,j)
            Pz(i+7*NB,j+13*NB)=sqrt(1./3)*p0m(i,j)
! 8v8v
            sum1=(0.,0.)
            sum2=(0.,0.)
            sum3=(0.,0.)
            sum4=(0.,0.)
            DO l=1,NB
               sum1=sum1
     &             +.5*(kz(i,l)*((2*g1m(l,j)-g2m(l,j)-3*g3m(l,j))/2)
     &                 +((2*g1m(i,l)-g2m(i,l)-3*g3m(i,l))/2)*kz(l,j))
               sum2=sum2+.5*(kz(i,l)*g2m(l,j)+g2m(i,l)*kz(l,j))
               sum3=sum3+.5*(kz(i,l)*g3m(l,j)+g3m(i,l)*kz(l,j))
               sum4=sum4
     &             +.5*(kz(i,l)*((2*g1m(l,j)+g2m(l,j)+3*g3m(l,j))/2)
     &                 +((2*g1m(i,l)+g2m(i,l)+3*g3m(i,l))/2)*kz(l,j))
            ENDDO
            Px(i+8*NB,j+8*NB)=-2*beta*(g1m(i,j)+g2m(i,j))*kx
            Px(i+8*NB,j+9*NB)=2*sqrt(3.)*beta*sum3
            Px(i+8*NB,j+10*NB)=2*sqrt(3.)*beta*g2m(i,j)*kx
     &                        -2*sqrt(3.)*beta*z*g3m(i,j)*ky
     &                        +z*ckm(i,j)
            Px(i+9*NB,j+9*NB)=-2*beta*(g1m(i,j)-g2m(i,j))*kx
            Px(i+9*NB,j+11*NB)=2*sqrt(3.)*beta*g2m(i,j)*kx
     &                        -2*sqrt(3.)*beta*z*g3m(i,j)*ky
     &                        -z*ckm(i,j)
            Px(i+10*NB,j+10*NB)=Px(i+9*NB,j+9*NB)
            Px(i+10*NB,j+11*NB)=-Px(i+8*NB,j+9*NB)
            Px(i+11*NB,j+11*NB)=Px(i+8*NB,j+8*NB)
            Py(i+8*NB,j+8*NB)=
     &           -beta*(2*g1m(i,j)-g2m(i,j)+3*g3m(i,j))*ky
     &           -sqrt(3.)/4*ckm(i,j)
            Py(i+8*NB,j+9*NB)=-2*sqrt(3.)*beta*z*sum2
            Py(i+8*NB,j+10*NB)=-sqrt(3.)*beta*(g2m(i,j)+g3m(i,j))*ky
     &                        -2*sqrt(3.)*beta*z*g3m(i,j)*kx
     &                        +.25*ckm(i,j)
            Py(i+9*NB,j+9*NB)=
     &           -beta*(2*g1m(i,j)+g2m(i,j)-3*g3m(i,j))*ky
     &           +3*sqrt(3.)/4*ckm(i,j)
            Py(i+9*NB,j+11*NB)=-sqrt(3.)*beta*(g2m(i,j)+g3m(i,j))*ky
     &                        -2*sqrt(3.)*beta*z*g3m(i,j)*kx
     &                        -.25*ckm(i,j)
            Py(i+10*NB,j+10*NB)=
     &           -beta*(2*g1m(i,j)+g2m(i,j)-3*g3m(i,j))*ky
     &           -3*sqrt(3.)/4*ckm(i,j)
            Py(i+10*NB,j+11*NB)=-Py(i+8*NB,j+9*NB)
            Py(i+11*NB,j+11*NB)=
     &           -beta*(2*g1m(i,j)-g2m(i,j)+3*g3m(i,j))*ky
     &           +sqrt(3.)/4*ckm(i,j)
            Pz(i+8*NB,j+8*NB)=-2*beta*sum1
            Pz(i+8*NB,j+9*NB)=
     &             2*sqrt(3.)*beta*(g3m(i,j)*kx-z*g2m(i,j)*ky)
     &            -.25*z*ckm(i,j)
            Pz(i+8*NB,j+10*NB)=-sqrt(3.)*beta*(sum2-sum3)
            Pz(i+8*NB,j+11*NB)=-3*sqrt(3.)/4*z*ckm(i,j)
            Pz(i+9*NB,j+9*NB)=-2*beta*sum4
            Pz(i+9*NB,j+10*NB)=sqrt(3.)/4*z*ckm(i,j)
            Pz(i+9*NB,j+11*NB)=Pz(i+8*NB,j+10*NB)
            Pz(i+10*NB,j+10*NB)=Pz(i+9*NB,j+9*NB)
            Pz(i+10*NB,j+11*NB)=
     &            -2*sqrt(3.)*beta*(g3m(i,j)*kx-z*g2m(i,j)*ky)
     &            -.25*z*ckm(i,j)
            Pz(i+11*NB,j+11*NB)=Pz(i+8*NB,j+8*NB)
! 8v7v
            Px(i+8*NB,j+12*NB)=-sqrt(6.)*beta*sum3
            Px(i+8*NB,j+13*NB)=-2*sqrt(6.)*beta*g2m(i,j)*kx
     &                         +2*sqrt(6.)*beta*z*g3m(i,j)*ky
     &                         -1./(2*sqrt(2.))*z*ckm(i,j)
            Px(i+9*NB,j+12*NB)=-2*sqrt(2.)*beta*g2m(i,j)*kx
     &                        +.5*sqrt(3./2)*z*ckm(i,j)
            Px(i+9*NB,j+13*NB)=3*sqrt(2.)*beta*sum3
            Px(i+10*NB,j+12*NB)=Px(i+9*NB,j+13*NB)
            Px(i+10*NB,j+13*NB)=-Px(i+9*NB,j+12*NB)
            Px(i+11*NB,j+12*NB)=2*sqrt(6.)*beta*g2m(i,j)*kx
     &                         +2*sqrt(6.)*beta*z*g3m(i,j)*ky
     &                         +1./(2*sqrt(2.))*z*ckm(i,j)
            Px(i+11*NB,j+13*NB)=Px(i+8*NB,j+12*NB)
            Py(i+8*NB,j+12*NB)=sqrt(6.)*beta*z*sum2
            Py(i+8*NB,j+13*NB)=
     &                2*sqrt(3./2)*beta*(g2m(i,j)+g3m(i,j))*ky
     &               +2*sqrt(6.)*beta*z*g3m(i,j)*kx
     &               +sqrt(1./2)*ckm(i,j)
            Py(i+9*NB,j+12*NB)=
     &                2*sqrt(1./2)*beta*(g2m(i,j)-3*g3m(i,j))*ky
            Py(i+9*NB,j+13*NB)=-3*sqrt(2.)*beta*z*sum2
            Py(i+10*NB,j+12*NB)=-Py(i+9*NB,j+13*NB)
            Py(i+10*NB,j+13*NB)=-Py(i+9*NB,j+12*NB)
            Py(i+11*NB,j+12*NB)=
     &               -2*sqrt(3./2)*beta*(g2m(i,j)+g3m(i,j))*ky
     &               +2*sqrt(6.)*beta*z*g3m(i,j)*kx
     &               +sqrt(1./2)*ckm(i,j)
            Py(i+11*NB,j+13*NB)=-Py(i+8*NB,j+12*NB)
            Pz(i+8*NB,j+12*NB)=
     &               -sqrt(6.)*beta*(g3m(i,j)*kx-z*g2m(i,j)*ky)
     &               -1./(2*sqrt(2.))*z*ckm(i,j)
            Pz(i+8*NB,j+13*NB)=2*sqrt(3./2)*beta*(sum2-sum3)
            Pz(i+9*NB,j+12*NB)=2*sqrt(1./2)*beta*(sum2+3*sum3)
            Pz(i+9*NB,j+13*NB)=
     &                3*sqrt(2.)*beta*(g3m(i,j)*kx-z*g2m(i,j)*ky)
     &                +.5*sqrt(3./2)*z*ckm(i,j)
            Pz(i+10*NB,j+12*NB)=
     &                3*sqrt(2.)*beta*(g3m(i,j)*kx+z*g2m(i,j)*ky)
     &                +.5*sqrt(3./2)*z*ckm(i,j)
            Pz(i+10*NB,j+13*NB)=-Pz(i+9*NB,j+12*NB)
            Pz(i+11*NB,j+12*NB)=-Pz(i+8*NB,j+13*NB)
            Pz(i+11*NB,j+13*NB)=
     &               -sqrt(6.)*beta*(g3m(i,j)*kx+z*g2m(i,j)*ky)
     &               -1./(2*sqrt(2.))*z*ckm(i,j)
! 7v7v
            Px(i+12*NB,j+12*NB)=-2*beta*g1m(i,j)*kx
            Px(i+13*NB,j+13*NB)=-2*beta*g1m(i,j)*kx
            Py(i+12*NB,j+12*NB)=-2*beta*g1m(i,j)*ky
            Py(i+13*NB,j+13*NB)=-2*beta*g1m(i,j)*ky
            sum1=(0.,0.)
            DO l=1,NB
               sum1=sum1+.5*(kz(i,l)*g1m(l,j)+g1m(i,l)*kz(l,j))
            ENDDO
            Pz(i+12*NB,j+12*NB)=-2*beta*sum1
            Pz(i+13*NB,j+13*NB)=-2*beta*sum1
         ENDDO
      ENDDO
      DO i=1,14*NB-1
         DO j=i+1,14*NB
            Px(j,i)=conjg(Px(i,j))
            Py(j,i)=conjg(Py(i,j))
            Pz(j,i)=conjg(Pz(i,j))
         ENDDO
      ENDDO
      END!end p110
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BME(NB,kz,pot8c,pot7c,pot6c,pot8v,pot7v,
     &               dmm,p0m,p1m,qm,gcm,g1m,g2m,g3m,ckm)
      PARAMETER (nzx=300)
      REAL L0,L1,L2,L3,sz(nzx),pot8c(NB,NB),pot7c(NB,NB),pot6c(NB,NB),
     &     pot8v(NB,NB),pot7v(NB,NB),dmm(NB,NB),p0m(NB,NB),p1m(NB,NB),
     &     qm(NB,NB),gcm(NB,NB),g1m(NB,NB),g2m(NB,NB),g3m(NB,NB),
     &     ckm(NB,NB)
      COMPLEX z,kz(NB,NB)
      COMMON/block02/e0,e1,d0,d1,dm,p0,p1,q,
     &               gLc,gL1,gL2,gL3,gc,g1,g2,g3,ck
      COMMON/block03/e0l,e1l,d0l,d1l,dml,p0l,p1l,ql,
     &               gLcl,gL1l,gL2l,gL3l,gcl,g1l,g2l,g3l,ckl
      COMMON/block04/e0r,e1r,d0r,d1r,dmr,p0r,p1r,qr,
     &               gLcr,gL1r,gL2r,gL3r,gcr,g1r,g2r,g3r,ckr
      COMMON/block05/L0,L1,L2,L3
      COMMON/block06/nz,sz,dz
      COMMON/block10/phi,eE0,tdu,eEz
      z=cmplx(0.,1.)
      vboffset=0.35
      DO i=1,NB
         DO j=i,NB
            sum1=0.
            DO l=1,nz
               x=sz(l)
               sum1=sum1+dz*basis(i,x)*dbasis(j,x)
            ENDDO
            kz(i,j)=-z*sum1
            kz(j,i)=conjg(kz(i,j))
         ENDDO
      ENDDO
      DO i=1,NB
         DO j=i,NB
            sum1=0.
            sum2=0.
            sum3=0.
            sum4=0.
            sum5=0.
            DO l=1,nz
               x=sz(l)
               IF (x.lt.(L1-L0/2)) THEN
                  sum1=sum1+dz*basis(i,x)*basis(j,x)
     &                     *(e0+(e0l-e0)*(1-vboffset)+(e1l-e0l)+d1l)
     &			   +dz*basis(i,x)*basis(j,x)*eEz*x
                  sum2=sum2+dz*basis(i,x)*basis(j,x)
     &                     *(e0+(e0l-e0-eEz*x)*(1-vboffset)+(e1l-e0l))
     &			 +dz*basis(i,x)*basis(j,x)*eEz*x
                  sum3=sum3+dz*basis(i,x)*basis(j,x)
     &                        *(e0+(e0l-e0-eEz*x)*(1-vboffset))
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
                  sum4=sum4+dz*basis(i,x)*basis(j,x)
     &                        *(-(e0l-e0)*vboffset)
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
                  sum5=sum5+dz*basis(i,x)*basis(j,x)
     &                        *(-(e0l-e0-eEz*x)*vboffset-d0l)
     &			 +dz*basis(i,x)*basis(j,x)*eEz*x
               ELSE IF (x.gt.(L1+L2-L0/2)) THEN
                  sum1=sum1+dz*basis(i,x)*basis(j,x)
     &                     *(e0+(e0r-e0)*(1-vboffset)+(e1r-e0r)+d1r)
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
                  sum2=sum2+dz*basis(i,x)*basis(j,x)
     &                     *(e0+(e0r-e0)*(1-vboffset)+(e1r-e0r))
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
                  sum3=sum3+dz*basis(i,x)*basis(j,x)
     &                        *(e0+(e0r-e0)*(1-vboffset))
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
                  sum4=sum4+dz*basis(i,x)*basis(j,x)
     &                        *(-(e0r-e0)*vboffset)
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
                  sum5=sum5+dz*basis(i,x)*basis(j,x)
     &                        *(-(e0r-e0)*vboffset-d0r)
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
               ELSE
                  sum1=sum1+dz*basis(i,x)*basis(j,x)*(e1+d1)
     &			 +dz*basis(i,x)*basis(j,x)*eEz*x
                  sum2=sum2+dz*basis(i,x)*basis(j,x)*e1
     &			 +dz*basis(i,x)*basis(j,x)*eEz*x
                  sum3=sum3+dz*basis(i,x)*basis(j,x)*e0
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
                  sum4=sum4+dz*basis(i,x)*basis(j,x)*0.
     &			 +dz*basis(i,x)*basis(j,x)*eEz*x
                  sum5=sum5+dz*basis(i,x)*basis(j,x)*(-d0)
     &			+dz*basis(i,x)*basis(j,x)*eEz*x
               ENDIF
            ENDDO
            pot8c(i,j)=sum1
            pot8c(j,i)=pot8c(i,j)
            pot7c(i,j)=sum2
            pot7c(j,i)=pot7c(i,j)
            pot6c(i,j)=sum3
            pot6c(j,i)=pot6c(i,j)
            pot8v(i,j)=sum4
            pot8v(j,i)=pot8v(i,j)
            pot7v(i,j)=sum5
            pot7v(j,i)=pot7v(i,j)
         ENDDO
      ENDDO
      DO i=1,NB
         DO j=i,NB
            sum1=0.
            sum2=0.
            sum3=0.
            sum4=0.
            sum5=0.
            sum6=0.
            sum7=0.
            sum8=0.
            sum9=0.
            DO l=1,nz
               x=sz(l)
               IF (x.lt.(L1-L0/2)) THEN
                  sum1=sum1+dz*basis(i,x)*basis(j,x)*dml
                  sum2=sum2+dz*basis(i,x)*basis(j,x)*p0l
                  sum3=sum3+dz*basis(i,x)*basis(j,x)*p1l
                  sum4=sum4+dz*basis(i,x)*basis(j,x)*ql
                  sum5=sum5+dz*basis(i,x)*basis(j,x)*gcl
                  sum6=sum6+dz*basis(i,x)*basis(j,x)*g1l
                  sum7=sum7+dz*basis(i,x)*basis(j,x)*g2l
                  sum8=sum8+dz*basis(i,x)*basis(j,x)*g3l
                  sum9=sum9+dz*basis(i,x)*basis(j,x)*ckl
               ELSE IF (x.gt.(L1+L2-L0/2)) THEN
                  sum1=sum1+dz*basis(i,x)*basis(j,x)*dmr
                  sum2=sum2+dz*basis(i,x)*basis(j,x)*p0r
                  sum3=sum3+dz*basis(i,x)*basis(j,x)*p1r
                  sum4=sum4+dz*basis(i,x)*basis(j,x)*qr
                  sum5=sum5+dz*basis(i,x)*basis(j,x)*gcr
                  sum6=sum6+dz*basis(i,x)*basis(j,x)*g1r
                  sum7=sum7+dz*basis(i,x)*basis(j,x)*g2r
                  sum8=sum8+dz*basis(i,x)*basis(j,x)*g3r
                  sum9=sum9+dz*basis(i,x)*basis(j,x)*ckr
               ELSE
                  sum1=sum1+dz*basis(i,x)*basis(j,x)*dm
                  sum2=sum2+dz*basis(i,x)*basis(j,x)*p0
                  sum3=sum3+dz*basis(i,x)*basis(j,x)*p1
                  sum4=sum4+dz*basis(i,x)*basis(j,x)*q
                  sum5=sum5+dz*basis(i,x)*basis(j,x)*gc
                  sum6=sum6+dz*basis(i,x)*basis(j,x)*g1
                  sum7=sum7+dz*basis(i,x)*basis(j,x)*g2
                  sum8=sum8+dz*basis(i,x)*basis(j,x)*g3
                  sum9=sum9+dz*basis(i,x)*basis(j,x)*ck
               ENDIF
            ENDDO
            dmm(i,j)=sum1
            dmm(j,i)=dmm(i,j)
            p0m(i,j)=sum2
            p0m(j,i)=p0m(i,j)
            p1m(i,j)=sum3
            p1m(j,i)=p1m(i,j)
            qm(i,j)=sum4
            qm(j,i)=qm(i,j)
            gcm(i,j)=sum5
            gcm(j,i)=gcm(i,j)
            g1m(i,j)=sum6
            g1m(j,i)=g1m(i,j)
            g2m(i,j)=sum7
            g2m(j,i)=g2m(i,j)
            g3m(i,j)=sum8
            g3m(j,i)=g3m(i,j)
            ckm(i,j)=sum9
            ckm(j,i)=ckm(i,j)
         ENDDO
      ENDDO
      END!end bme
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION basis(n,x)
      PARAMETER (pi=3.14159)
      REAL L0,L1,L2,L3
      COMMON/block05/L0,L1,L2,L3
      basis=sqrt(2./L0)*sin(n*pi/L0*(x+L0/2))
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION dbasis(n,x)
      PARAMETER (pi=3.14159)
      REAL L0,L1,L2,L3
      COMMON/block05/L0,L1,L2,L3
      dbasis=sqrt(2./L0)*(n*pi/L0)*cos(n*pi/L0*(x+L0/2))
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE param_AlGaAs(x,tem,e0,e1,d0,d1,dm,p0,p1,q,
     &                        gLc,gL1,gL2,gL3,gc,g1,g2,g3,ck)
      COMMON/block01/hb,xm0,beta
      e0=1.519+1.087*x+0.438*x**2     !eV
      e0=e0-5.405e-4*tem**2/(tem+204)      !temperature dependence
      e1=4.488+0.052*x       !eV
      d0=0.341-0.041*x       !eV
      d1=0.171-0.021*x       !eV
      dm=-0.05               !eV
      p0=1.0493-0.1523*x     !eV.nm
      p1=0.478               !eV.nm
      q=0.8165               !eV.nm
      gLc=1/(0.0665+0.0835*x)
      gL1=6.85-3.6*x
      gL2=2.1-1.45*x
      gL3=2.9-1.69*x
      ck=-.00034+0.00054*x   !eV.nm
      gc=gLc-(2*xm0)/(3*hb*hb)*(
     &         p0*p0*(2/e0+1/(e0+d0))-p1*p1*(1/(e1-e0)+2/(e1+d1-e0))
     &        +(4./3)*p0*p1*dm*(1/(e0*(e1+d1-e0))-1/((e0+d0)*(e1-e0))))
      g1=gL1-(2*xm0)/(3*hb*hb)*(p0*p0/e0+q*q/(e1+d1)+q*q/e1
     &                         +(2./3)*p0*p1*dm/(e0*(e1+d1)))
      g2=gL2-(2*xm0)/(3*hb*hb)*(.5*(p0*p0/e0-q*q/e1)
     &                         +(1./3)*p0*p1*dm/(e0*(e1+d1)))
      g3=gL3-(2*xm0)/(3*hb*hb)*(.5*(p0*p0/e0+q*q/e1)
     &                         +(1./3)*p0*p1*dm/(e0*(e1+d1)))
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INCLUDE 'chiev.f'
