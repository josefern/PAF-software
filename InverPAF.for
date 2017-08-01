c**********************************************************************
c                                                                     *
c                Inver P A F    (25/03/2017)                          * 
c                                                                     *
c     A software tool for inversion of terrain deformations           *
c        as free-geometry extended bodies for anomalous pressure.     *
c                                                                     *
c     A.G. Camacho (1), J. Fernández (1), F. Cannavó (2)              *
c     (1) Institute of Geosciences (CSIC-UCM), Madrid, Spain          *
c     (2) Osservatorio Etneo, INGV, Catania, Italy.                   *
c                                                                     *
c     See PAFmanual.txt for description of the input and output files *
c      and operation approach                                         *
c                                                                     *
c**********************************************************************
      use ifport
      use ifqwin

      implicit real*8(a-h,o-z)
      parameter (ms=5000,mc=350000,mb=2000,mm=5000)
      type (xycoord) xy
      TYPE (rccoord) curpos

      dimension fs(2),x(ms),y(ms),z(ms),rab(ms),
     -dg(ms),wg(ms),rg(ms),
     -dz(ms),wz(ms),zcp(ms),zcd(ms),rz(ms),zo(ms),za(ms),
     -dx(ms),wx(ms),xcp(ms),xcd(ms),rx(ms),xo(ms),xa(ms),
     -dy(ms),wy(ms),ycp(ms),ycd(ms),ry(ms),yo(ms),ya(ms),
     -qg(mc),qz(mc),xc(mc),yc(mc),zc(mc),vc(mc), zas(ms),di(ms),sz(ms),
     -cdz(ms),cdx(ms),cdy(ms),cpz(ms),cpx(ms),cpy(ms),cdg(ms),cpg(ms),
     -jb(mb),xb(mb),yb(mb),sen(mc)
      integer*4 m(mc),p(mc),nor(mc),js(2),ipa(mc),ifa(mc)
      integer*2 nada,idz(ms),idh(ms),ks(2),ld(50),lp(50)
      character*50 obs,mod,gri,fil,flos,tex(12),raya
      character*9 hoy,hora1,hora2, fic(mm)
      character*130 texto

      data raya/'-'/,pi/3.1141592/
      call seed(2016)
      idum=1000
      blun=9.0
      nuf=1
      sig=0.25
      dmu=10
      depli=3000.         

c---------------------------Input data and covariance matrix-------------!
      j=time()
      idum=j-int(j/1.e3)*1.e3
      iu=1
      iu0=0
      ih0=0
      los=0
      dpr=0
      is3=0
      lx=560
      ly=190
      lz=460
      tx=360
      ty=360
      lx1=lx-tx/2
      lx2=lx+tx/2
      lx3=lx-tx/2+380
      lx4=lx+tx/2+380
      ly1=ly-ty/2
      ly2=ly+ty/2

      open(1,file='CellsConfig.txt')
      read(1,'(///a7)') hora2
      read(1,*) pem
      read(1,*) dlap
      read(1,*) seli
      read(1,*) idib
      read(1,*) nepo
      read(1,*) iep1
      read(1,'(////a9)') hora2
      if(nc.gt.mc) call dim(mc)
                                          
      iz=1
      ih=1

      open(2,file='DeforData.txt')
      i=0
      ax=9.d9
      ay=ax
      bx=-ax
      by=-ay
      xm=0
      ym=0
      zm=0
    1 read(2,*,err=2,end=2) xp,yp,zp
      if(xp.eq.0.and.yp.eq.0) go to 2
      if(xp.lt.ax) ax=xp
      if(xp.gt.bx) bx=xp
      if(yp.lt.ay) ay=yp
      if(yp.gt.by) by=yp
      xm=xm+xp
      ym=ym+yp
      zm=zm+zp
      i=i+1
      if(i.gt.ms) write(*,*) ' ***** Num of data site > ',i6
      if(i.gt.ms) stop
      x(i)=xp
      y(i)=yp
      z(i)=zp
      go to 1
    2 ns=i
      close(2)
      xm=xm/ns
      ym=ym/ns
      zm=zm/ns
      if(ns.eq.0) stop
      na=4
      if(ns.gt.500) na=9
      if(ns.lt.30) na=2
      fe=(bx-ax+by-ay)/2./30000.
      do 3 i=1,ns
      x(i)=x(i)-xm
      y(i)=y(i)-ym
    3 z(i)=z(i)-zm
                                          
      do 6 i=1,nepo
      no=iep1-1+i
      if(no.le.99999) write(fic(i),'(i5)') no
      if(no.le.9999) write(fic(i),'(a1,i4)') '0',no
      if(no.le.999) write(fic(i),'(a2,i3)') '00',no
      if(no.le.99) write(fic(i),'(a3,i2)') '000',no
      if(no.le.9) write(fic(i),'(a4,i1)') '0000',no
    6 continue
                                          
      tic=4000*fe
      i=log10(tic)
      tic=nint(tic/10**i)*10**i
      zme=0
      siz=0
      j=0
    5 read(1,*,err=27,end=27) xp,yp,zp,px,py,pz
      if(xp.eq.0.and.yp.eq.0) go to 27
      j=j+1
      p(j)=0
      sen(j)=0
      zc(j)=zp-zm                        
      zme=zme+zc(j)
      xc(j)=xp-xm
      yc(j)=yp-ym
      vo=px*py*pz/fe/fe/fe
      r=(px+py+pz)/3./fe
      siz=siz+r
      if(r.lt.500) r=500
      r2=r*r
      at2=0
      dv=0.
      sw=9.d9
         cc=0
      do 4 i=1,ns
      cx=(x(i)-xc(j))/fe                         
      cy=(y(i)-yc(j))/fe
      cz=(z(i)-zc(j))/fe
      d2=cx*cx+cy*cy+cz*cz
      if(d2.lt.sw) dv=cz
      if(d2.lt.sw) sw=d2
         if(d2.lt.r2) d2=r2
      c=cz*cz/d2/d2/d2
      at2=at2+c
         if(c.gt.cc) cc=c
    4 continue
         sen(j)=sqrt((at2-cc)/(ns-1))*1.e9
         if(dv.lt.depli) sen(j)=sen(j) *dv/depli
      at2=at2/ns
      qz(j)=vo*vo*at2                    
      vc(j)=vo*fe*fe*fe                  
      go to 5
   27 close(1)
      nc=j
      zme=zme/nc
      siz=siz/nc

      if(idib.eq.1) then
      open(1,file='map.bln')                              
      k=1                                                 
    9 read(1,*,end=8) j                                   
      do 7 i=k,(j+k-1)                                    
      read(1,*) xb(i),yb(i)                               
      xb(i)= (xb(i)-xm)/fe/120+lx                         
      yb(i)=-(yb(i)-ym)/fe/120+ly                         
      jb(i)=0                                             
      if(i.ge.mb) go to 8                                 
    7 continue                                            
      jb(k)=1                                             
      k=k+j                                               
      go to 9                                             
    8 nb=k-1                                              
      close(1)                                            
      endif


c-------------Inversion process--------------------------------------

      do 88 ii=1,nepo                                      
      write(obs,'(a1,a5,a4)') 'D',fic(ii),'.txt'           
      write(mod,'(a1,a5,a4)') 'M',fic(ii),'.txt'           

      if(nepo.eq.1) write(obs,'(a13)') 'DeforData.txt'    
      if(nepo.eq.1) write(mod,'(a10)') 'ModPAF.txt'       
           
      nz=0
      nx=0
      ny=0
      swz=0
      swx=0
      swy=0
      sw=0
      ez=1
      ex=1
      ey=1
      if(iz.eq.0) ez=0
      if(ih.eq.0) ex=0
      if(ih.eq.0) ey=0
      open(2,file=obs)
      k=0
      do 13 i=1,ns
      read(2,*,err=14,end=14) xp,yp,zp,(tex(j),j=1,6)
      k=k+1
      if(xp.eq.0.and.yp.eq.0.) go to 14
      idz(i)=0
      dz(i)=0.
      if(iz.ne.1.or.tex(1).eq.raya.or.tex(2).eq.raya) go to 12
       read(tex(1),*) dz(i)
       read(tex(2),*) erz
       if(erz.eq.0) go to 12
        idz(i)=1
        wz(i)= erz
        swz=swz+wz(i)
        nz=nz+1
   12 idh(i)=0
      dx(i)=0.
      dy(i)=0.
      if(ih.eq.0) go to 13
      if(tex(3).eq.raya.or.tex(4).eq.raya.or.tex(5).eq.raya) go to 13
       read(tex(3),*) dx(i)
       read(tex(4),*) erx
       read(tex(5),*) dy(i)
       read(tex(6),*) ery
       if(ery.eq.0.and.erx.eq.0) go to 13
       idh(i)=1
       wx(i)= erx   !1
       wy(i)= ery   !1
       swx=swx+wx(i)
       swy=swy+wy(i)
       nx=nx+1
       ny=ny+1
   13 continue
   14 close(2)
      if(k.ne.ns) write(*,*) 'Data error'
      if(k.ne.ns) stop
      if(nz.gt.0) swz=swz/nz
      if(nx.gt.0) swx=swx/nx
      if(ny.gt.0) swy=swy/ny
      cx=swx
      cy=swy
      cz=swz

      sdz=0
      sdx=0
      sdy=0
      emp=0.
      swu=0
      swz=0
      swx=0
      swy=0
      dzm=0
      do 15 i=1,ns
      dzm=dzm+dz(i)
      if(idz(i).eq.1) wz(i)=wz(i)/cz
      if(idz(i).eq.2) wz(i)=wz(i)/cz2
      if(idz(i).eq.3) wz(i)=wz(i)/cz3
      if(cx.ne.0) wx(i)=wx(i)/cx
      if(cy.ne.0) wy(i)=wy(i)/cy
      if(sw.ne.0) wg(i)=wg(i)/sw
      zcp(i)=0.
      zcd(i)=0.
      xcp(i)=0.
      xcd(i)=0.
      ycp(i)=0.
      ycd(i)=0.
      if(idz(i).eq.1) then
       wz(i)=ez*wz(i)
       sdz=sdz+dz(i)*dz(i)
       emp=emp+dz(i)*dz(i)*wz(i)
       swu=swu+wz(i)
       swz=swz+wz(i)
      endif
       dx(i)=dx(i)-dx0
       dy(i)=dy(i)-dy0
       if(idh(i).eq.1) then
       wx(i)=ex*wx(i)
       wy(i)=ey*wy(i)
       sdx=sdx+dx(i)*dx(i)
       sdy=sdy+dy(i)*dy(i)
       emp=emp+dx(i)*dx(i)*wx(i)
       swu=swu+wx(i)
       swx=swx+wx(i)
       emp=emp+dy(i)*dy(i)*wy(i)
       swu=swu+wy(i)
       swy=swy+wy(i)
      endif
   15 continue
      dzm=dzm/ns
      if(nz.gt.0) sdz=sqrt(sdz/nz)
      if(ii.eq.1) tdef=nint(sdz*2.)/1.0
      if(nx.gt.0) sdx=sqrt(sdx/nx)
      if(ny.gt.0) sdy=sqrt(sdy/ny)
      if(swu.gt.0.) emp1=sqrt(emp/swu)
      nrej=0 
      c=5 
      do 16 i=1,ns
      cx=0
      cy=0
      cz=0
      s=0
      do 17 j=1,ns    
      if(i.eq.j) go to 17
      d=1./((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
      s=s+d
      cx=cx+dx(j)*d
      cy=cy+dy(j)*d
      cz=cz+dz(j)*d
   17 continue   
      cx=cx/s    
      cy=cy/s    
      cz=cz/s    
      if(abs(dz(i)-cz).gt.c*sdz) idz(i)=0
      if(abs(dx(i)-cx).gt.c*sdx) idh(i)=0
      if(abs(dy(i)-cy).gt.c*sdy) idh(i)=0
      k=idz(i)*idh(i)
   16 if(k.eq.0) nrej=nrej+1
      
      if(idib.eq.1) then
c------------------------------------------------------------Graphics--!
      nada=settextcolorrgb(#000000)                                    !
      nada=setbkcolorrgb(#f0f0f0)                                      !
      call clearscreen($clearscreen)                                   !
      nada=initializefonts()                                           !
      nada=setcolorrgb(#000000)                                        !
      nada=setfont('t''Times New Roman''h20w9')                        !
      write(texto,'(a)') 'PAF-Inversion of site displacements'         !
      call moveto(780,10,xy)                                           !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'as 3D free-geometry pressurized bodies'      !
      call moveto(800,35,xy)                                           !
      nada=setfont('t''Times New Roman''h16w7')                        !      
      call outgtext(texto)                                             !
      nada=setfont('t''Times New Roman''h12w5')                        !
      write(texto,'(a)') 'A.G.Camacho, 2016'                           !
      k=lx+tx/2.+4                                                     !
      call moveto(k,120,xy)                                            !
      call setgtextrotation(900)                                       !
      call outgtext(texto)                                             !
      call setgtextrotation(000)                                       !
      call date(hoy)                                                   !
      call time(hora1)                                                 !
      nada=setcolorrgb(#000000)                                        !
      nada=setfont('t''Courier''h16w8')                                !
      write(texto,'(a)') hoy                                           !
      call moveto(900,65,xy)                                           !
      call outgtext(texto)                                             !
      write(texto,'(2a)') 'Start ..... ',hora1                         !
      call moveto(850,90,xy)                                           ! 
      call outgtext(texto)                                             !
      call ventana(lx1,lz-50,lx2,lz+100)                               !  
      call ventana(lx3,lz-50,lx4,lz+100)                               !  
      call ventana(20,lz-50,330,lz+100)                                !  
      call ventana(lx1,ly1,lx2,ly2)                                    !  
      nada=setcolorrgb(#dfffdf)                                        ! 
      nada=rectangle($gfillinterior,lx1,lz,lx2,lz+100)                 ! 
      nada=rectangle($gfillinterior,lx3,lz,lx4,lz+100)                 ! 
      nada=rectangle($gfillinterior,20,lz,330,lz+100)                  ! 
      nada=rectangle($gfillinterior,lx1,ly1,lx2,ly2)                   ! 
      nada=setcolorrgb(#00a000)                                        !
      call ventana( 65,70,110,85)                                      !
      call ventana(130,70,250,85)                                      !
      nada=setcolorrgb(#000000)                                        !
      write(texto,'(a,a)') 'Data file:  ',obs                          !
      call moveto(60,25,xy)                                            !
      call outgtext(texto)                                             !
      write(texto,'(a,5x,a)') 'Number','Scatt.(cm)'                    !
      call moveto(70,50,xy)                                            !
      call outgtext(texto)                                             !
      write(texto,'(i6,2x,3f5.1)') ns,sdz,sdx,sdy                      !
      call moveto(60,70,xy)                                            !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Up   EW   SN '                               !
      call moveto(145,88,xy)                                           !
      call outgtext(texto)                                             !
      
      do 19 j=1,nc                                                     !
      i=p(j)                                                           !
      if(i.le.0) go to 19                                              !
      kx= xc(j)/fe/120+lx                                              !
      ky=-yc(j)/fe/120+ly                                              !
      kz=-zc(j)/fe/120+lz                                              !
      nada=setcolorrgb(#90f090)                                        !
        if(sen(j).lt.seli) nada=setcolorrgb(#c0ffc0)                   !
        kk=vc(j)**0.33/120/3.                                          !
      if(kx.lt.lx1.or.kx.gt.lx2.or.ky.gt.ly2.or.ky.lt.ly1) go to 19    !
        nada=rectangle($gfillinterior,kx-kk,ky-kk,kx+kk,ky+kk)         !
      if((kz+kk).lt.(lz+100))                                          !  
     -  nada=rectangle($gfillinterior,kx-kk,kz-kk,kx+kk,kz+kk)         !
   19 continue                                                         !
      nada=setcolorrgb(#00a000)                                        !
      kk=1                                                             !
      do 10 i=1,nb                                                     !
      kx=xb(i)                                                         !
      ky=yb(i)                                                         !
      if(kk.eq.1) jb(i)=1                                              !
      if(kk.eq.0) kk=1                                                 !
      if(kx.lt.lx1.or.kx.gt.lx2.or.ky.gt.ly2.or.ky.lt.ly1) go to 10    !
      if(jb(i).eq.1) call moveto(kx,ky,xy)                             !
      if(jb(i).eq.0) nada=lineto(kx,ky)                                !
      kk=0                                                             !
   10 continue                                                         !
      nada=setcolorrgb(#000000)                                        !
      c=tic/120/fe                                                     !
      do 11 i=1,12                                                     !
      k=(i-1)*tic/120/fe                                               ! 
      call moveto(lx1,lz+k,xy)                                         !
      if((lz+k).lt.(lz+100)) nada=lineto(lx1+4,lz+k)                   !
      call moveto(lx2,lz+k,xy)                                         !
      if((lz+k).lt.(lz+100)) nada=lineto(lx2-4,lz+k)                   !
      call moveto(lx3,lz+k,xy)                                         !
      if((lz+k).lt.(lz+100)) nada=lineto(lx3+4,lz+k)                   ! 
      call moveto(lx4,lz+k,xy)                                         !
      if((lz+k).lt.(lz+100)) nada=lineto(lx4-4,lz+k)                   !
      call moveto(lx1,ly1+k,xy)                                        !
      if((ly1+k).lt.ly2) nada=lineto(lx1+4,ly1+k)                      !
      call moveto(lx2,ly1+k,xy)                                        !
      if((ly1+k).lt.ly2) nada=lineto(lx2-4,ly1+k)                      !
      call moveto(lx1+k,ly2,xy)                                        !
      if((lx1+k).lt.lx2) nada=lineto(lx1+k,ly2-4)                      !
      call moveto(lx1+k,ly1,xy)                                        !
      if((lx1+k).lt.lx2) nada=lineto(lx1+k,ly1+4)                      !
      iy=(i-1)*c                                                       !
      if(iy.gt.110) go to 11                                           !
      call moveto(18,lz+iy,xy)                                         !
      nada=lineto(22,lz+iy)                                            !
   11 continue                                                         !
      call moveto(16,lz-50+iy,xy)                                      !
      nada=lineto(20,lz-50+iy)                                         !
      do 18 i=1,ns                                                     !
      kx= x(i)/fe/120+lx                                               !
      ky=-y(i)/fe/120+ly                                               !
      kz=-z(i)/fe/120+lz                                               !
      nada=setcolorrgb(#000000)                                        !
      if(kx.le.lx1.or.kx.ge.lx2.or.ky.ge.ly2.or.ky.le.ly1) go to 18    !      
      nada=ellipse($gborder,kx-2,ky-2,kx+2,ky+2)                       !
      nada=ellipse($gborder,kx-2,kz-2,kx+2,kz+2)                       !
      kx= y(i)/fe/120+lx+380                                           !
      nada=ellipse($gborder,kx-2,kz-2,kx+2,kz+2)                       !
   18 continue                                                         !
      nada=setcolorrgb(#000000)                                        !
      write(texto,'(a,a)')                                             !
     -' Press. Mod.   RMS Resid.(cm)   Depth'                          !
      call moveto(20,225,xy)                                           !
      call outgtext(texto)                                             !
      write(texto,'(a,a)')                                             !
     -' Step  (MPa)   Up   WE   SN     (m)'                            !
      call moveto(25,241,xy)                                           !
      call outgtext(texto)                                             !
      call ventana(16,262,320,279)                                     !
      write(texto,'(a)') '- Fit'                                       !
      call moveto(150,300,xy)                                          !
      nada=setcolorrgb(#009000)                                        !
      call outgtext(texto)                                             !
      call ventana(lx1,lz+140,lx2,lz+280)                              !
      nada=setcolorrgb(#000000)                                        !
      write(texto,'(10x,a)') ' Initial steps'                          !
      call moveto(40,264,xy)                                           !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Model growth -->'                            !
      call moveto(20,lz+108,xy)                                        !
      call outgtext(texto)                                             !
      iy=390                                                           !
      ix=780                                                           !
      nada=setcolorrgb(#dfffdf)                                        !
      nada=rectangle($gfillinterior,ix+100,225,ix+225,277)             !
      nada=rectangle($gfillinterior,ix+190,150,ix+205,165)             !
      nada=setcolorrgb(#000000)                                        ! 
      write(texto,'(a)') 'Survey points..'                             !
      call moveto(ix+70,150,xy)                                        !
      call outgtext(texto)                                             !
      nada=ellipse($gborder,ix+195,155,ix+200,160)                     !
      write(texto,'(a)') 'Press. cells '                               !
      call moveto(ix+90,190,xy)                                        !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Signific..'                                  !
      call moveto(ix+10,230,xy)                                        !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'No signif.'                                  !
      call moveto(ix+10,255,xy)                                        !
      call outgtext(texto)                                             !
      write(texto,'(a)') ' -     +   Prev.'                            !
      call moveto(ix+105,212,xy)                                       !
      call outgtext(texto)                                             !
      write(texto,'(a,i4)')   'Signif.limit ....',nint(seli)           !
      call moveto(ix+70,300,xy)                                        !
      call outgtext(texto)                                             !
      call colorfu(1,1,1,ix+115,235,ix+120,240)                        !
      nada=setcolorrgb(#fff090)                                        !
      nada=rectangle($gfillinterior,ix+115,260,ix+120,265)             !
      nada=setcolorrgb(#90f090)                                        !
      nada=rectangle($gfillinterior,ix+205,235,ix+210,240)             !
      nada=setcolorrgb(#c0ffc0)                                        !
      nada=rectangle($gfillinterior,ix+205,260,ix+210,265)             !
      call colorfu(1,2,1,ix+160,235,ix+165,240)                        !
      nada=setcolorrgb(#60f0f0)                                        !
      nada=rectangle($gfillinterior,ix+160,260,ix+165,265)             !
      nada=setcolorrgb(#00b0f0)                                        !
      nada=rectangle($gfillinterior,65-1,620,65,645)                   !
      nada=setcolorrgb(#f0a000)                                        !
      nada=rectangle($gfillinterior,144-1,630,144,645)                 !
      nada=setcolorrgb(#0000f0)                                        !
      nada=rectangle($gfillinterior,230-1,610,230,645)                 !
      nada=setcolorrgb(#a0f0f0)                                        !
      nada=rectangle($gfillinterior,320-1,610,320,645)                 !
      nada=setcolorrgb(#000000)                                        !
      write(texto,'(a)') 'Obs..     Mod..    Out ..    no dy ..'       !
      call moveto(20,628,xy)                                           !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Model - Plan view'                           !
      call moveto(lx-70,377,xy)                                        !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Model - WE elevation view'                   !
      call moveto(lx-70,568,xy)                                        !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Model - SN elevation view'                   !
      call moveto(lx+310,568,xy)                                       !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Defor.Fit - WE elevation view'               !
      call moveto(lx-90,752,xy)                                        !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Defor.Fit - SN elevation view'               !
      call moveto(lx+290,752,xy)                                       !
      call outgtext(texto)                                             !
      call setgtextrotation(900)                                       !
      write(texto,'(a,f4.1,a)') 'Tic=',tdef,' cm'                      !
      call moveto(358,710,xy)                                          !
      call outgtext(texto)                                             !
      write(texto,'(a,i6,a)') 'Tic=',nint(tic),' m'                    !
      call moveto(358,240,xy)                                          !
      call outgtext(texto)                                             !
      call moveto(358,560,xy)                                          !
      call outgtext(texto)                                             !
      call moveto(  5,560,xy)                                          !
      call outgtext(texto)                                             !
      write(texto,'(a)') 'Press scal'                                  !
      call moveto(6,390,xy)                                            !
      call outgtext(texto)                                             !
      call setgtextrotation(000)                                       !
      write(texto,'(a,i5)') 'Mean cell size(m)..',nint(siz)            !
      call moveto(60,130,xy)                                           !
      call outgtext(texto)                                             !
      write(texto,'(a,f5.1)') 'Press.contr.(MPa)..',pem                !
      call moveto(60,150,xy)                                           !
      call outgtext(texto)                                             !
      write(texto,'(a,i5)')   'Smooth.Coeff.......',nint(dlap)         !
      call moveto(60,170,xy)                                           !
      call outgtext(texto)                                             !
 
      call ventana(lx1,lz+140,lx2,lz+280)                              !
      call ventana(lx3,lz+140,lx4,lz+280)                              !      
      nada=setcolorrgb(#000000)                                        !
      kk=lz+210+dzm*30./tdef                                           !
      call moveto(lx1,kk,xy)                                           !
      nada=lineto(lx2,kk)                                              ! 
      call moveto(lx3,kk,xy)                                           !
      nada=lineto(lx4,kk)                                              ! 
      do 25 i=1,10                                                     !
      j=kk+30*i                                                        !
      call moveto(lx1,j,xy)                                            !
      if(j.lt.lz+280) nada=lineto(lx1+4,j)                             !
      call moveto(lx3,j,xy)                                            !
      if(j.lt.lz+280) nada=lineto(lx3+4,j)                             !
      j=kk-30*i                                                        !
      call moveto(lx1,j,xy)                                            !
      if(j.gt.lz+140) nada=lineto(lx1+4,j)                             !
      call moveto(lx3,j,xy)                                            !
      if(j.gt.lz+140) nada=lineto(lx3+4,j)                             !
   25 continue                                                         !
      do 26 i=1,ns                                                     ! 
      ix= x(i)/fe/120+lx                                               !
      nada=setcolorrgb(#00b0f0)                                        !
      if(idz(i).eq.0) nada=setcolorrgb(#0000f0)                        !
      if(idh(i).eq.0) nada=setcolorrgb(#0000f0)                        !
      call moveto(ix,kk,xy)                                            !
      k= ix+dx(i)*30./tdef                                             !
      iy=kk-dz(i)*30./tdef                                             !
      nada=lineto(k,iy)                                                !
      ix= y(i)/fe/120+lx +380                                          !
      call moveto(ix,kk,xy)                                            !
      k= ix+dy(i)*30./tdef                                             !
      iy=kk-dz(i)*30./tdef                                             !
      if(dy(i).eq.0.) nada=setcolorrgb(#a0f0f0)                        !
      nada=lineto(k,iy)                                                !
   26 continue                                                         !
c----------------------------------------------------------------------!
      endif
      nu5=2
      nrem=0
      svp=0.
      fps1=0
      ozs=0
      oxs=0
      oys=0
      ol2=0
      ol3=0
      nu=0
      npr=0
      do 28 j=1,nc
      p(j)=-j
   28 nor(j)=0
          nsig=0
      kpf=1
      ls=1
   20 continue
      fitu=999.d9
      ncs=nc
      if(na.gt.1.and.nu.gt.10) ncs=nc/na
      do 21 i=1,ns
      if(idz(i).eq.1) then
        zo(i)=dz(i) -ozs
        zas(i)=zcp(i)
      endif
      if(idz(i).eq.2) then
        zo(i)=dz(i)-ol2
        zas(i)=-xcp(i)*cs1+ycp(i)*ss1+zcp(i)*cl1
      endif
      if(idz(i).eq.3) then
        zo(i)=dz(i) -ol3
        zas(i)=-xcp(i)*cs2+ycp(i)*ss2+zcp(i)*cl2
      endif
      if(idh(i).eq.0) go to 21
        xo(i)=dx(i) -oxs
        yo(i)=dy(i) -oys
   21 continue
              dlp=dlap*1.e-4/siz *(1.-1./(nu+1))     
      do 30 ja=1,ncs
      j=ja
      if(ncs.lt.nc) call RANDOM_NUMBER(ranval)
      if(ncs.lt.nc) j=(nc-1)*ranval+1
      vo=vc(j)
      sv=svp+qz(j)*dlp
      if(p(j).ge.0) go to 30
      call deformod(dmu,sig,fra,xc(j),yc(j),zc(j),vo,ms,ns,x,y,z,1,1,
     -cdz,cdx,cdy,cdg,cpz,cpx,cpy,cpg)
      do 24 k=1,2
      pb=2*k-3
      std=0.
      sdd=0.
      do 23 i=1,ns
      xp=cpx(i)
      yp=cpy(i)
      zp=cpz(i)
      if(idz(i).ne.0) then
       if(idz(i).eq.1) c=zp
       if(idz(i).eq.2) c=-xp*cs1+yp*ss1+zp*cl1
       if(idz(i).eq.3) c=-xp*cs2+yp*ss2+zp*cl2
       za(i)=zas(i)+c*pb
       zp=za(i)
       ww=wz(i)*zp
       std=std+ww*zo(i)
       sdd=sdd+ww*zp
      endif
      if(idh(i).eq.0) go to 23
       xa(i)=xcp(i)+xp*pb
       xp=xa(i)
       ww=wx(i)*xp
       std=std+ww*xo(i)
       sdd=sdd+ww*xp
       ya(i)=ycp(i)+yp*pb
       yp=ya(i)
       ww=wy(i)*yp
       std=std+ww*yo(i)
       sdd=sdd+ww*yp
   23 continue
      fp=std/(sdd+sv)
      if(fp.le.0.) go to 24
      em=0.
      do 22 i=1,ns
      if(idz(i).eq.0) go to 22
      c=zo(i)-za(i)*fp
      em=em+c*c*wz(i)
      if(idh(i).eq.0) go to 22
      c=xo(i)-xa(i)*fp
      em=em+c*c*wx(i)
      c=yo(i)-ya(i)*fp
      em=em+c*c*wy(i)
   22 continue
      em=em/swu
      c=em+fp*fp*sv
      if(c.ge.fitu) go to 24
       fitu=c
       fs(1)=fp
       js(1)=j
       ks(1)=k
       kps=1
       emp=em
   24 continue
   30 continue
      j=js(ls)
      if(j.eq.0) go to 50
      vo=vc(j)

       call deformod(dmu,sig,fra,xc(j),yc(j),zc(j),vo,ms,ns,x,y,z,
     - 1,1,cdz,cdx,cdy,cdg,cpz,cpx,cpy,cpg)
       fitu1=fitu
       npr=npr+1
       emp1=emp
       fps1=fs(1)
       jp=js(1)
       p(jp)=ks(1)
       nor(jp)=fps1*10
       if(npr.le.10) lp(npr)=jp
       pb=2*ks(1)-3
       svp=svp+qz(jp)*dlp
       ipa(jp)=1
       ifa(jp)=1
       sdx=0.
       sdy=0.
       sdz=0.
       sdl2=0
       sdl3=0
       cz2=0.
       cz3=0.
       cz=0
       cx=0
       cy=0

       do 33 i=1,ns
       xp=cpx(i)
       yp=cpy(i)
       zp=cpz(i)
       gp=cpg(i)
       zcp(i)=zcp(i)+zp*pb
       xcp(i)=xcp(i)+xp*pb
       ycp(i)=ycp(i)+yp*pb
       d=zcd(i)
       if(idz(i).eq.2) d=-xcd(i)*cs1+ycd(i)*ss1+zcd(i)*cl1
       if(idz(i).eq.3) d=-xcd(i)*cs2+ycd(i)*ss2+zcd(i)*cl2
       zz=ozs
       if(idz(i).eq.2) zz=ol2
       if(idz(i).eq.3) zz=ol3
       c=dz(i)-zz
       d=zcp(i)
       if(idz(i).eq.2) d=-xcp(i)*cs1+ycp(i)*ss1+zcp(i)*cl1
       if(idz(i).eq.3) d=-xcp(i)*cs2+ycp(i)*ss2+zcp(i)*cl2
       rz(i)=c-d*fps1
       c=c-d
       if(idz(i).eq.1) then
         sdz=sdz+wz(i)*rz(i)*rz(i)
         cz=cz+wz(i)*rz(i)
       endif
       if(idz(i).eq.2) then
         sdl2=sdl2+wz(i)*rz(i)*rz(i)
         cz2=cz2+wz(i)*rz(i)
       endif
       if(idz(i).eq.3) then
         sdl3=sdl3+wz(i)*rz(i)*rz(i)
         cz3=cz3+wz(i)*rz(i)
       endif
       rx(i)=dx(i)-oxs -xcp(i)*fps1
       ry(i)=dy(i)-oys -ycp(i)*fps1
       if(idh(i).eq.0) go to 33
       sdx=sdx+wx(i)*rx(i)*rx(i)
       sdy=sdy+wy(i)*ry(i)*ry(i)
       cx=cx+wx(i)*rx(i)
       cy=cy+wy(i)*ry(i)
   33  continue
       if(swz.gt.0.) sdz=sqrt(sdz/swz)
       if(swl2.gt.0.) sdl2=sqrt(sdl2/swl2)
       if(swl3.gt.0.) sdl3=sqrt(sdl3/swl3)
       if(swx.gt.0.) sdx=sqrt(sdx/swx)
       if(swy.gt.0.) sdy=sqrt(sdy/swy)
       if(swz2.gt.0) cz2=cz2/swl2
       if(swz3.gt.0) cz3=cz3/swl3
       if(swx.gt.0.) cx=cx/swx
       if(swy.gt.0.) cy=cy/swy
       if(swz.gt.0.) cz=cz/swz
       if(iu0.eq.1) ozs=ozs+cz*0.5
       if(ih0.eq.1) oxs=oxs+cx*0.5
       if(ih0.eq.1) oys=oys+cy*0.5
       if(il0.eq.1) then
        ol2=ol2+cz2*0.5
        ol3=ol3+cz3*0.5
       endif
      
      tt=fps1/pem
      call dmedian(ms,ns,rz,em)
      em=em/0.6745*blun *tt*tt*tt
      mout=0
      d=0.
      do 34 i=1,ns
   34 if(wz(i).gt.d) d=wz(i)
      do 35 i=1,ns
      c=rz(i)/em
      if(c.le.1.) go to 35
       c=c*c-1.
       if(c.lt.0.5) mout=mout+1
       c=d/(1.+c*c)
       if(c.lt.wz(i)) wz(i)=c
   35 continue
           if(sen(jp).ge.seli) nsig=nsig+1
      nu=npr

      if(idib.eq.0) go to 36
c-----------------------------------------------------Graphics---------!
      write(texto,'(i5,f7.1,1x,3f5.1,f9.0)')                           !
     -npr,fps1,sdz,sdx,sdy,zc(jp)+zm                                   !
      nada=setcolorrgb(#ffffff)                                        !
      nada=rectangle($gfillinterior,18,262,320,279)                    !
      nada=setcolorrgb(#000000)                                        !
      call moveto(20,264,xy)                                           !
      call outgtext(texto)                                             !
      ix=20+nu*2                                                       !
      kx= xc(jp)/fe/120+lx                                             !
      ky=-yc(jp)/fe/120+ly                                             !
      kz=-zc(jp)/fe/120+lz                                             !
         c=sqrt(fps1)                                                  !
      if(nu.eq.1) fep=c                                                !
      if(nu.eq.1) fiu=fitu                                             !
      if(ks(ls).eq.1) nada=setcolorrgb(#f00000)                        !
      if(ks(ls).eq.2) nada=setcolorrgb(#0000f0)                        !
                if(sen(jp).lt.seli) then                               !
                    if(ks(ls).eq.1) nada=setcolorrgb(#fff090)          !
                    if(ks(ls).eq.2) nada=setcolorrgb(#60f0f0)          !
                endif                                                  !
      k=390-c/fep*150                                                  !
      if(k.lt.290) k=290                                               !      
      if(nu.le.154) nada=rectangle($gfillinterior,ix,390,ix+2,k)       ! 
      kk=vc(jp)**0.33/fe/120./1.9                                      ! 
      iy=lz-zc(js(ls))/fe/120.                                         !
      if(nu.lt.154) nada=rectangle($gfillinterior,ix,iy-kk,ix+2,iy+kk) !
      if(kx.gt.lx1.and.kx.lt.lx2.and.ky.lt.ly2.and.ky.gt.ly1) then     !
        nada=rectangle($gfillinterior,kx-kk,ky-kk,kx+kk,ky+kk)         !
        i=lz+100                                                       !
        if((kz+kk).lt.i) i=kz+kk                                       ! 
        if((kz-kk).lt.(lz+100))                                        !
     -   nada=rectangle($gfillinterior,kx-kk,kz-kk,kx+kk,i)            !
        k=yc(jp)/fe/120+lx+380                                         !
        if((kz-kk).lt.(lz+100))                                        !
     -   nada=rectangle($gfillinterior,k-kk,kz-kk,k+kk,i)              !
      endif                                                            !
      k=390-fitu/fiu*30                                                !
      if(k.lt.290) k=290                                               !      
      nada=setcolorrgb(#009000)                                        !
      if(nu.le.154) nada=rectangle($gfillinterior,ix,k,ix+1,k+1)       !
      nu5=nu5+1                                                        !        
      if(nu5.lt.10) go to 36                                           !
      nu5=0                                                            !      
      call ventana(lx1,lz+140,lx2,lz+280)                              !
      call ventana(lx3,lz+140,lx4,lz+280)                              !      
      nada=setcolorrgb(#000000)                                        !
      kk=lz+210+dzm*30./tdef                                           !
      call moveto(lx1,kk,xy)                                           !
      nada=lineto(lx2,kk)                                              ! 
      call moveto(lx3,kk,xy)                                           !
      nada=lineto(lx4,kk)                                              ! 
      do 38 i=1,10                                                     !
      j=kk+30*i                                                        !
      call moveto(lx1,j,xy)                                            !
      if(j.lt.lz+280) nada=lineto(lx1+4,j)                             !
      call moveto(lx3,j,xy)                                            !
      if(j.lt.lz+280) nada=lineto(lx3+4,j)                             !
      j=kk-30*i                                                        !
      call moveto(lx1,j,xy)                                            !
      if(j.gt.lz+140) nada=lineto(lx1+4,j)                             !
      call moveto(lx3,j,xy)                                            !
      if(j.gt.lz+140) nada=lineto(lx3+4,j)                             !
   38 continue                                                         !
      do 37 i=1,ns                                                     !
      nada=setcolorrgb(#00b0f0)                                        !
      if(idz(i).eq.0.or.idh(i).eq.0) nada=setcolorrgb(#0000f0)         !
      ix= x(i)/fe/120+lx                                               !
      iy= y(i)/fe/120+lx+380                                           !
      call moveto(ix,kk,xy)                                            !
      j= ix+dx(i)*30./tdef                                             !
      k=kk-dz(i)*30./tdef                                              !
      nada=lineto(j,k)                                                 !
      call moveto(iy,kk,xy)                                            !
      l= iy+dy(i)*30./tdef                                             !
      if(dy(i).eq.0.) nada=setcolorrgb(#a0f0f0)                        !
      nada=lineto(l,k)                                                 !
      if(idz(i).eq.0.or.idh(i).eq.0) go to 37                          !
      nada=setcolorrgb(#f0a000)                                        !
      call moveto(ix,kk,xy)                                            ! 
      j=j-rx(i)*30./tdef                                               !
      k=k+rz(i)*30./tdef                                               !      
      if(k.ge.lz+280.or.k.le.lz+140) go to 37                          !  
      nada=lineto(j,k)                                                 !
      call moveto(iy,kk,xy)                                            ! 
      l=l-ry(i)*30./tdef                                               !
      nada=lineto(l,k)                                                 !
   37 continue                                                         !
c-----------------------------------------------------------------------
   36 continue
      l=ls
      ls=l
      js(1)=0
      js(2)=0
      if(fps1.gt.0.and.fps1.lt.pem.and.npr.gt.2) kpf=0
      if(kpf.ne.0) go to 20

c----------------------------------------------------------------------+
c--------------  End of process ---------------------------------------+

   50 continue
      call time(hora2)
      if(idib.eq.1) then
c---------------------------------------------------------Graphics------
      nada=setcolorrgb(#000000)                                        !
      write(texto,'(a,a)') 'End ....... ',hora2                        !
      call moveto(850,110,xy)                                          !
      call outgtext(texto)                                             !
      call ventana(990,350,1035,365)                                   !
      nada=setcolorrgb(#000000)                                        !
      write(texto,'(a,i7)') 'Num.sign.cells  ',nsig                    !
      call moveto(850,350,xy)                                          !
      call outgtext(texto)                                             !
      call ventana(220,680,265,695)                                    !
      call ventana(220,715,265,730)                                    !
      nada=setcolorrgb(#000000)                                        !
      write(texto,'(a,f7.3)') 'Misfit value ...',fitu1                 !
      call moveto(80,680,xy)                                           !
      call outgtext(texto)                                             !
      write(texto,'(a,i6)') 'Reject. data ...',nrej                    !
      call moveto(80,715,xy)                                           !
      call outgtext(texto)                                             !
      if(nu.eq.0) then                                                 !
      write(texto,*) '*** Error: error in the data. No cells. '        !
      call moveto(320,260,xy)                                          !
      call outgtext(texto)                                             !
      endif                                                            !
      write(texto,'(a)') ' ===> Press (Enter)'                         !
      call moveto(100,760,xy)                                          !
      call outgtext(texto)                                             !
      read(*,'(i1)') i                                                 !
c-----------------------------------------------------------------------
      endif

      open(1,file=mod)
      open(2,file='CellsConfig.txt')
      do 59 i=1,12
      read(2,'(a130)') texto
   59 write(1,'(a130)') texto
      read(2,'(a130)') texto
      read(2,'(a130)') texto
      read(2,'(a130)') texto
      write(1,'(/a//a/a)')
     -' Cells: location (UTM, m) , sides (m)  and  press (MPa)',
     -'    X      Y       Z      sx   sy   sz     Press   Signi',
     -'---------------------------------------------------------'
      nn=0
      do 51 j=1,nc
   51 if(p(j).gt.0.or.m(j).gt.0) nn=nn+1

      totmas=0
      totpre=0
      swz=0
        npp=0
        npn=0
      do 52 j=1,nc
      read(2,*) jx,jy,jz,kx,ky,kz
      vol=kx*ky*kz
      pb=0.
      i=p(j)
      if(i.le.0) go to 52
      pb=(2*i-3)*fps1    !presion
        if(i.eq.1) npn=npn+1
        if(i.eq.2) npp=npp+1
   58 totpre=totpre+vol*abs(pb)
      swz=swz+jz*vol*abs(pb)
      k=sen(j)
      write(1,'(i7,i8,i7,2x,3i5,f9.2,i5)') jx,jy,jz,kx,ky,kz,pb,k
   52 continue
      swz=swz/totpre
      write(1,'(a,a9,a)') 
     -'--------------------- Date:',hoy,'------------------'
      write(1,202) ns,nc,npn+npp,npn,npp,sig,dmu,na
  202 format(/' Num. data points =',i6
     -/' Num. total cells=',i6
     -/' Num.filled cells=',i5,4x,'Neg, Pos =',2i5
     -/' Medium param.:',6x,'Poisson=',f4.2,4x,'Share Mod.=',f4.0,'GPa'
     -/' Random explor.coeff.=',i4)
      write(1,203) fps1,totpre/1.e9,swz,sdz,sdx,sdy,fitu1
  203 format(' Pressure model:',6x,'Press. contrast (-+)=',f7.2,' MPa'
     -/22x,'Press*vol:',f7.1,' MPA*Km3',4x,' Mean mod.depth =',f7.0
     -/' RMS residuals (cm):  Up=',f5.2,4x,'WE=',f5.2,4x,'SN=',f5.2,
     -5x,'Misfit=',f9.4)
   53 close(2)
      write(1,'(a,2a12)') ' Initial and final exec.times:',hora1,hora2

      write(1,'(///20x,a)') 'Observed, modeled, and residual values'
      write(1,'(/a,10x,a,12x,a,16x,a/2x,a,6x,a,7x,a,8x,a,10x,a,10x,a
     -/a,a)')
     -'Data point loc(UTM, m)','dz (cm)','dx(cm)','dy(cm)','X','Y',
     -'Z','obs  mod  res','obs  mod  res','obs mod res    Weight'
     -,'----------------------------------------------',
     -'-----------------------------------------------'
      k=0
      emg=0
      emu=0
      do 55 i=1,ns
      if(idz(i).eq.0) go to 55
      emu=emu+dz(i)*dz(i)
      k=k+1
   55 continue
      if(k.gt.0) emu=sqrt(emu/k)
      cg=0
      cz=0
      c2=0
      c3=0
      do 56 i=1,ns
      xp=xm+x(i)
      yp=ym+y(i)
      zp=zm+z(i)
      ie=0
      if(emg.ne.0.) ie=abs(rg(i))/emg
      if(d.lt.0.01) d=0.01
      up=0
      if(emu.ne.0.) up=abs(rz(i))/emu
      zcal=dz(i)-rz(i)
      gcal=dg(i)-rg(i)
      xcal=dx(i)-rx(i)
      ycal=dy(i)-ry(i)
      d=0
      if(idz(i).eq.2) d=dl2+dxl2*x(i)/1000.+dyl2*y(i)/1000.
      if(idz(i).eq.3) d=dl3+dxl3*x(i)/1000.+dyl3*y(i)/1000.
      write(1,210) nint(xp),nint(yp),nint(zp),
     -dz(i)+d,zcal+d,rz(i),
     -dx(i)+dx0,xcal+dx0,rx(i),dy(i)+dy0,ycal+dy0,ry(i),wz(i)
  210 format(i6,i8,i5,2x,3(3f7.2,1x),f5.1)
   56 continue
      write(1,'(a,a//60x,a,i4)') 
     -'----------------------------------------------',
     -'-----------------------------------------------',
     -'Number of outlier values=',mout
      close(1)

   88 continue                                  !<<<

c      close(6)                                 
      stop
      end

c***********************************************************************
c   Dimension                                                          *
c***********************************************************************
      subroutine dim(m)
      write(*,200) m
  200 format(/6x,'*** Error: data size > max. dimension',i6,' !!!'/
     -10x,'Reduce the data side or change dimension in the code'/)
      stop
      end
C***********************************************************************
c Surface modeled elevation changes (cm)                               *
c  for pressure and mass = 1 MPa, 1 Kg                                 *
c***********************************************************************
      subroutine deformod(dmu,sig,fra,xc,yc,zc,vo,ms,ns,x,y,z,kt,kf,
     -cdz,cdx,cdy,cdg,cpz,cpx,cpy,cpg)
      implicit real*8(a-h,o-z)
      dimension x(ms),y(ms),z(ms),
     -cdz(ms),cdx(ms),cdy(ms),cpz(ms),cpx(ms),cpy(ms),cdg(ms),cpg(ms),
     -cx(14),cy(14),cz(14)
      data cx/0, 1, 0, 0, 0.7071, 0.7071,0.7071,-0.7071, 0.0000, 0.0000,
     - 0.57735,-0.57735, 0.57735, 0.57735/
      data cy/0, 0, 1, 0, 0.7071,-0.7071,0.0000, 0.0000, 0.7071,-0.7071,
     - 0.57735, 0.57735,-0.57735, 0.57735/
      data cz/0, 0, 0, 1, 0.0000, 0.0000,0.7071, 0.7071, 0.7071, 0.7071,
     - 0.57735, 0.57735, 0.57735,-0.57735/
      csp=0.577
      pi=3.14159265359
      dera=pi/180.
      u1=9.8066/4./pi/dmu/1.e8          
      ra3=0.750/pi*vo
      siz=vo**0.333333
      u2=(1.-sig)/dmu*ra3/10.   
      do 10 i=1,ns
      xp=x(i)-xc
      yp=y(i)-yc
      zp=z(i)-zc
      x2=xp*xp
      y2=yp*yp
      z2=zp*zp
      d2=x2+y2+z2 +siz
      d=sqrt(d2)
      d3=d2*d
      cdx(i)=0
      cdy(i)=0
      cdz(i)=0
      cdg(i)=0
      cpx(i)=0
      cpy(i)=0
      cpz(i)=0
      cpg(i)=0
      if(kt.eq.0) then   
        cdz(i)=-u1/d*(2.*(1.-sig)+z2/d2)*vo        
        pp=-u1/d*(zp/d2-(1.-2.*sig)/(d+zp))*vo     
        cdx(i)=pp*xp                               
        cdy(i)=pp*yp                               
        cdg(i)=vo*zp/d3*6.672d-3 - cdz(i)*fra      
      endif
      if(kf.le.1) go to 2
      if(kt.eq.2.and.kf.gt.1) then
        if(kf.eq.2.and.yp.gt.0) go to 10
        if(kf.eq.3.and.yp.lt.0) go to 10
        if(kf.eq.4.and.zp.gt.0) go to 10
        if(kf.eq.5.and.zp.lt.0) go to 10
      endif
      if(kt.eq.3.and.kf.gt.1) then
        if(kf.eq.2.and.xp.gt.0) go to 10
        if(kf.eq.3.and.xp.lt.0) go to 10
        if(kf.eq.4.and.zp.gt.0) go to 10
        if(kf.eq.5.and.zp.lt.0) go to 10
      endif
      if(kt.eq.4.and.kf.gt.1) then
        if(kf.eq.2.and.xp.gt.0) go to 10
        if(kf.eq.3.and.xp.lt.0) go to 10
        if(kf.eq.4.and.yp.gt.0) go to 10
        if(kf.eq.5.and.yp.lt.0) go to 10
      endif
      if(kt.eq.5.and.kf.gt.1) then
        s=xp-yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.zp.gt.0) go to 10
        if(kf.eq.5.and.zp.lt.0) go to 10
      endif
      if(kt.eq.6.and.kf.gt.1) then
        s=xp+yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.zp.gt.0) go to 10
        if(kf.eq.5.and.zp.lt.0) go to 10
      endif
      if(kt.eq.7.and.kf.gt.1) then
        s=xp-zp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.yp.gt.0) go to 10
        if(kf.eq.5.and.yp.lt.0) go to 10
      endif
      if(kt.eq.8.and.kf.gt.1) then
        s=xp+zp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.yp.gt.0) go to 10
        if(kf.eq.5.and.yp.lt.0) go to 10
      endif
      if(kt.eq.9.and.kf.gt.1) then
        s=yp-zp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.xp.gt.0) go to 10
        if(kf.eq.5.and.xp.lt.0) go to 10
      endif
      if(kt.eq.10.and.kf.gt.1) then
        s=yp+zp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.xp.gt.0) go to 10
        if(kf.eq.5.and.xp.lt.0) go to 10
      endif
      if(kt.eq.11.and.kf.gt.1) then        
        s=xp-yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        s=yp-zp
        if(kf.eq.4.and.s.gt.0) go to 10
        if(kf.eq.5.and.s.lt.0) go to 10
      endif
      if(kt.eq.12.and.kf.gt.1) then    
        s=xp+yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        s=yp-zp
        if(kf.eq.4.and.s.gt.0) go to 10
        if(kf.eq.5.and.s.lt.0) go to 10
      endif
      if(kt.eq.13.and.kf.gt.1) then     
        s=xp+yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        s=yp+zp
        if(kf.eq.4.and.s.gt.0) go to 10
        if(kf.eq.5.and.s.lt.0) go to 10
      endif
      if(kt.eq.14.and.kf.gt.1) then     
        s=xp-yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        s=yp+zp
        if(kf.eq.4.and.s.gt.0) go to 10
        if(kf.eq.5.and.s.lt.0) go to 10
      endif
    2 p=u2/d3  
      if(kt.eq.1) then    
       cpx(i)=p*xp *1.28
       cpy(i)=p*yp *1.28
       cpz(i)=p*zp
       go to 3
      endif
      if(kt.gt.1) p=p*1.0   
      if(kf.gt.1) p=p*1.0   
      s=p* d*(1-csp)        
      c=d2-cx(kt)*cx(kt)*x2-cy(kt)*cy(kt)*y2-cz(kt)*cz(kt)*z2
      c=p* (cx(kt)*xp+cy(kt)*yp+cz(kt)*zp)/sqrt(c)*csp
      cpx(i)=s*cx(kt)+c*(1-cx(kt))*xp
      cpy(i)=s*cy(kt)+c*(1-cy(kt))*yp
      cpz(i)=s*cz(kt)+c*(1-cz(kt))*zp
    3 cpg(i)=-cpz(i)*fra      
   10 continue
      return
      end

c****************************************************************
c     Square window                                             *
c****************************************************************
      subroutine ventana(i1,j1,i2,j2)
      use portlib
      use msflib
      nada=setcolorrgb(#909090)
      nada=rectangle($gfillinterior,i1-2,j1-2,i1-1,j2+1)
      nada=rectangle($gfillinterior,i1-2,j1-2,i2+1,j1-1)
      nada=setcolorrgb(#fefefe)
      nada=rectangle($gfillinterior,i1-1,j2+2,i2+2,j2+2)
      nada=rectangle($gfillinterior,i2+2,j1,i2+2,j2+2)
      nada=setcolorrgb(#ffffff)
      nada=rectangle($gfillinterior,i1,j1,i2,j2)
      return
      end
c************************************************************************
c  Color   b g r
c************************************************************************
      subroutine colorfu(i,ks,kf,ix1,iy1,ix2,iy2)
      use msflib
      integer*4 col(55)
      integer*2 kc(2,2,14)
      data kc/   1, 2,  1, 2,    
     -          29,30,  3, 4,    
     -          31,32,  5, 6,    
     -          33,34,  7, 8,    
     -          35,36,  9,10,    
     -          37,38, 11,12,    
     -          39,40, 13,14,    
     -          41,42, 15,16,    
     -          43,44, 17,18,    
     -          45,46, 19,20,    
     -          47,48, 21,22,    
     -          49,50, 23,24,    
     -          51,52, 25,26,    
     -          53,54, 27,28/    
      col( 1)=#f00000     
      col( 2)=#0000f0     
      col( 3)=#f00000
      col( 4)=#fff090
      col( 5)=#00b000
      col( 6)=#40c040
      col( 7)=#0000c0
      col( 8)=#5050c0
      col( 9)=#705000   
      col(10)=#908040   
      col(11)=#906020   
      col(12)=#708020   
      col(13)=#700060  
      col(14)=#a05080  
      col(15)=#702080  
      col(16)=#a02060  
      col(17)=#005060  
      col(18)=#408080  
      col(19)=#206080  
      col(20)=#208060  
      col(21)=#808080   
      col(22)=#707070   
      col(23)=#608080   
      col(24)=#806060   
      col(25)=#806080   
      col(26)=#608080   
      col(27)=#808060   
      col(28)=#606080   
      col(29)=#fff090   
      col(30)=#f0f080   
      col(31)=#00f0f0   
      col(32)=#60f0f0   
      col(33)=#f060f0   
      col(34)=#f0a0f0   
      col(35)=#70f0a0  
      col(36)=#a0f0b0  
      col(37)=#90f090  
      col(38)=#70f0b0  
      col(39)=#70a0a0  
      col(40)=#f0c0b0  
      col(41)=#f0a0b0  
      col(42)=#f0c0a0  
      col(43)=#70a0f0  
      col(44)=#a0c0f0  
      col(45)=#a0a0f0  
      col(46)=#70c0f0  
      col(47)=#b0b0b0 
      col(48)=#a0a0a0 
      col(49)=#a0b0b0 
      col(50)=#b0a0a0 
      col(51)=#b0a0b0 
      col(52)=#a0b0b0 
      col(53)=#b0b0a0 
      col(54)=#a0a0b0 
      jf=kf
      if(kf.gt.2) jf=2
      k=kc(ks,jf,i)
      nada=setcolorrgb(col(k))
      if(i.eq.0) nada=setcolorrgb(#b0b0b0)
      nada=rectangle($gfillinterior,ix1,iy1,ix2,iy2)
      return
      end
c***********************************************************************
      subroutine ORDEN(mx,nx,x,ix)
      real*8 x(mx)
      integer ix(mx)
      ix(1)=1
      if(nx.eq.1) return
      ia=1
      ib=nx
      do 1 i=1,nx
      if(x(i).gt.x(ib)) ib=i
    1 if(x(i).lt.x(ia)) ia=i
      ix(1)=ia
      do 2 i=2,nx
      i1=i-1
      im=ix(i1)
      ia=ib
      do 3 j=1,nx
      if(x(j).le.x(im)) go to 3
      if(x(j).lt.x(ia)) ia=j
    3 continue
      ix(i)=ia
    2 continue
      return
      end
c***********************************************************************
c   median of absolute values for n<m different numbers x(i)           c
c***********************************************************************
      subroutine dmedian(m,n,x,xa)
      implicit real*8(a-h,o-z)
      dimension x(m)
      xi=x(1)
      xs=xi
      do 1 i=1,n
      xx=x(i)
      if(xx.lt.xi)  xi=xx
      if(xx.gt.xs)  xs=xx
        do 6 j=i+1,n
    6   if(xx.eq.x(j)) x(j)=x(j)*1.00001
    1 continue
      j=0
    4 xa=xs
      do 2 i=1,n
      xx=x(i)
      if(xx.lt.xa.and.xx.gt.xi) xa=xx
    2 continue
      j=j+1
      if(j.eq.n) go to 9
      xi=xa
      xb=xi
      do 3 i=1,n
      xx=x(i)
      if(xx.gt.xb.and.xx.lt.xs) xb=xx
    3 continue
      j=j+1
      xa=(xa+xb)/2.
      xs=xb
      if(j.lt.n) go to 4
    9 return
      end
c***********************************************************************
c***************************** end code ********************************
