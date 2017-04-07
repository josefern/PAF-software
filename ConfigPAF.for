c**********************************************************************
c                                                                     *
c               ConfigPAF        (25/03/2017)                         *
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
      use ifqwin
      type (xycoord) xy
      parameter (ms=5000,mc=350000)
      real*4 at(mc),di(ms)
      integer x(ms),y(ms),z(ms),mada,
     -xb(mc),yb(mc),dx(mc),dy(mc),zt(mc),zb(mc)
      character*50 texto
      common ns,zm,fe,base,techo,pp,nc,pem,jsmo
      iex=1
      jsmo=200
      pem=0.5
      nada=settextcolorrgb(#000000)
      nada=setbkcolorrgb(#ffffff)
      nada=initializefonts()
      call clearscreen($clearscreen)
    1 open(1,file='DeforData.txt')
      xm=0.
      ym=0.
      zm=0.
      ax=9.d9
      bx=-ax
      ay=ax
      by=-ay
      i=1
    3 read(1,*,err=2,end=2) xx,yy,c
      if(xx.eq.0.) go to 2
      if(xx.lt.ax) ax=xx
      if(xx.gt.bx) bx=xx
      if(yy.lt.ay) ay=yy
      if(yy.gt.by) by=yy
      xm=xm+xx
      ym=ym+yy
      zm=zm+c
      i=i+1
      if(i.lt.ms) go to 3
      write(*,'(//20x,a,i4/)') ' Error: Num of stations > ',ms
      stop
    2 ns=i-1
      if(ns.gt.1) go to 4
      write(*,'(//20x,a,a12/)') ' Error in data file ',fobs
      stop
    4 rewind 1
      c=(bx-ax+by-ay)/2.
      fe=c/30000. 
      xm=nint(xm/ns/10)*10.
      ym=nint(ym/ns/10)*10.
      zm=nint(zm/ns/10)*10.
      techo=-9.d9     
      do 5 i=1,ns
      read(1,*) xx,yy,zz
      z(i)=(zz-zm)/fe
      if(z(i).gt.techo) techo=z(i)    
      x(i)=(xx-xm)/fe
    5 y(i)=(yy-ym)/fe
      close (1)
      ax=(ax-xm)/fe
      bx=(bx-xm)/fe
      ay=(ay-ym)/fe
      by=(by-ym)/fe
      pp=500
      base=-13000
      call di1()
      avm=pp*pp*pp*pp*1.d-10   
      nc=0
      kr=0
      nada=setfont('t''Courier''h17w7')
   10 ncr=nc
      kr=kr+1
      zs=techo
      ax=ax-17000
      ay=ay-17000
      bx=bx+17000
      by=by+17000
   12 nada=setcolorrgb(#ffffff)            
      nada=rectangle($gfillinterior,200,150,400,170)
      k=100-(zs-base)/(techo-base)*100
      call moveto(200,150,xy)
      nada=setcolorrgb(#000000)
      call outgtext(texto)
      tz=pp*0.8
      call especapa(zs,avm,ms,ns,x,y,z,tz) 
      if(tz.gt.1) go to 11
         zs=zs-pp*0.5
         if(zs.gt.base) go to 12
   11 if(kr.eq.2) tz=tz*2 
      if(kr.eq.3) tz=tz*2.7
      call celdascapa(ax,bx,ay,by,zs,tz,avm,ms,ns,x,y,z,
     -mc,ncr,xb,yb,zt,zb,dx,dy,at)
      zs=zs-tz                      
      if(zs.gt.base) go to 12  
      l=nc
      do 15 i=nc+1,ncr
      if(nc.eq.0) go to 16
      xx=dx(i) /2.       ! 1.9
      yy=dy(i) /2.       ! 1.9    
      pr=(zt(i)-zb(i)) /1.9
      k=0
      zz=(zt(i)+zb(i))/2.
      do 14 j=1,nc
      if(abs(xb(j)-xb(i)).gt.xx) go to 14
      if(abs(yb(j)-yb(i)).gt.yy) go to 14
      zj=(zt(j)+zb(j))/2.
      if(abs(zj-zz).gt.pr) go to 14
      k=1  
   14 continue
      if(k.eq.1) go to 15
   16 l=l+1
      xb(l)=xb(i)
      yb(l)=yb(i)
      zt(l)=zt(i)
      zb(l)=zb(i)
      dx(l)=dx(i)
      dy(l)=dy(i) 
   15 continue
      nc=l
      if(kr.lt.iex) go to 10   
   40 write(*,'(///a)') ' '
      call di2(ire)
      open(1,file='CellsConfig.txt')
      call step(ms,ns,x,y,di,c) 
      write(1,200)  pem,jsmo,5,1,1,0,' '
  200 format('      PAF-modeling of pressure sources'
     -//'   Parameters'
     -/' -----------------------------------------' 
     -/ f7.1,3x,'....Pressure change (MPa)        '
     -/   i7,3x,'....Smoothing coeff. (0<smo<1000)'
     -/   i7,3x,'....Significance limit (0<sig<10)'
     -/   i7,3x,'....Graphic output  0:no  1:YES  ' 
     -/   i7,3x,'....Number of data epochs        ' 
     -/   i7,3x,'....First epoch identification   ' 
     -/' -----------------------------------------' 
     -//2x,'Cells: location (UTM, m) and sides (m)',a)
      write(1,*) '    X      Y       Z      sx   sy   sz '   
      write(1,*) ' --------------------------------------'   
      do 41 i=1,nc
      ix=xb(i)*fe+xm
      iy=yb(i)*fe+ym
      iz=(zt(i)+zb(i))/2.*fe+zm  
      kx=dx(i)*fe
      ky=dy(i)*fe
      kz=(zt(i)-zb(i))*fe
   41 write(1,'(2i8,i7,2x,3i5,2x,3i5)') 
     -ix,iy,iz,kx,ky,kz   !,0,0,nint(at(i)*10)
      write(1,*) ' --------------------------------------'   
      write(1,'(5x,a,i7)') ' Number of cells=',nc
      close(1)
      stop
      end
c******************************************************************
c  Dialog window  1                                               *   
c******************************************************************
      subroutine di1()
      use dialogm
      include 'CellsConfig.fd'
      logical retlog,nada
      character*50 texto
      common ns,zm,fe,base,techo,pp,nc,pem,jsmo
      external fin
      type(dialog) dlg
      retlog=DlgInit(IDD_CellsConfig, dlg)
      retlog=DlgSet(dlg,IDC_TEXT_t1,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t2,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_BOX_b1 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e1,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e2,.false.,dlg_enable)
      retlog=DlgSetSub(dlg,IDCANCEL,fin)
      write(texto,'(i4)') ns
      nada=DlgSet(dlg,IDC_EDIT_e8,texto)
      write(texto,'(i7)') nint(base*fe+zm)
      nada=DlgSet(dlg,IDC_EDIT_e4,texto)
      write(texto,'(i7)') nint(techo*fe+zm)
      nada=DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(i6)') nint(pp*fe)
      nada=DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') jsmo
      nada=DlgSet(dlg,IDC_EDIT_e1,texto)
      write(texto,'(f3.1)') pem
      nada=DlgSet(dlg,IDC_EDIT_e2,texto)
      retint = DlgModal( dlg )
      nada=DlgGet(dlg,IDC_EDIT_e4,texto)
      read(texto,*) k
      base=(k-zm)/fe
      nada=DlgGet(dlg,IDC_EDIT_e6,texto)
      read(texto,*) k
      techo=(k-zm)/fe
      nada=DlgGet(dlg,IDC_EDIT_e5,texto)
      read(texto,*) k
      pp=k/fe
      return
      
      end
c******************************************************************
c  Dialog window 2                                                *   
c******************************************************************
      subroutine di2(ire)
      use dialogm
      include 'CellsConfig.fd'
      logical retlog,nada
      character*50 texto
      external fin
      common ns,zm,fe,base,techo,pp,nc,pem,jsmo
      type(dialog) dlg
              ! supress compiler warnings for unreferenced arguments
      if(.not.DlgInit(IDD_CellsConfig,dlg)) then
      write(*,*) "error: resource not found,  vale?"
      else
      retlog=DlgSet(dlg,IDC_TEXT_t4,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t5,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t9,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_EDIT_e4,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e5,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e6,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_BOX_b3 ,.false.,dlg_enable)
      write(texto,'(i3)') jsmo
      nada=DlgSet(dlg,IDC_EDIT_e1,texto)
      write(texto,'(f3.1)') pem
      nada=DlgSet(dlg,IDC_EDIT_e2,texto)
      write(texto,'(i4)') ns
      nada=DlgSet(dlg,IDC_EDIT_e8,texto)
      write(texto,'(f8.0)') base*fe+zm
      nada=DlgSet(dlg,IDC_EDIT_e4,texto)
      write(texto,'(f7.0)') techo*fe+zm
      nada=DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f6.0)') pp*fe
      nada=DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i6)') nc
      nada=DlgSet(dlg,IDC_EDIT_e9,texto)
      retlog=DlgSetSub(dlg,IDCANCEL,fin)
      retint = DlgModal( dlg )
      nada=DlgGet(dlg,IDC_EDIT_e1,texto)
      read(texto,*) jsmo
      nada=DlgGet(dlg,IDC_EDIT_e2,texto)
      read(texto,*) pem
      end if
      return
      end
c***********************************************************************
      subroutine fin
      use dialogm
      include 'CellsConfig.fd'
      stop
      end
c***********************************************************************
c   Median for n<m different numbers x(i)
c***********************************************************************
      subroutine dmedian(m,n,x,xa)
      real*4 x(m)
      xi=x(1)
      xs=x(1)
      do 1 i=1,n
      if(x(i).lt.xi)  xi=x(i)
      if(x(i).gt.xs)  xs=x(i)
    1 continue
      j=0
    4 xa=xs
      do 2 i=1,n
      if(x(i).lt.xa.and.x(i).gt.xi) xa=x(i)
    2 continue
      j=j+1
      if(j.eq.n) go to 9
      xi=xa
      xb=xi
      do 3 i=1,n
      if(x(i).gt.xb.and.x(i).lt.xs) xb=x(i)
    3 continue
      j=j+1
      xa=(xa+xb)/2.
      xs=xb
      if(j.lt.n) go to 4
    9 return
      end
c***********************************************************************
c     Sensitivity calculus
c***********************************************************************
      subroutine at2(ms,ns,x,y,z,xx,yy,zz,dx,dy,dz,se,av,zp)
      integer x(ms),y(ms),z(ms)
      av=0.
      t=(dx+dy+dz)/3
      if(t.lt.500) t=500   
      t2=t*t
      zp=0
      sp=0 
      do 1 i=1,ns
      zr=zz-z(i)
      xr=xx-x(i)
      yr=yy-y(i)
      d=xr*xr+yr*yr+zr*zr
      if(d.lt.t2) d=t2
      c=zr*zr/d/d/d
      av=av+c
      d=1/d/d
      sp=sp+d
      zp=zp+z(i)*d
    1 continue
      av=av/ns
      se=av*1.d16
      t2=dx*dy*dz
      av=av*t2*t2  
      zp=zp/sp
      return
      end
c***********************************************************************
c     Step dr for autocorrelation analysis
c***********************************************************************
      subroutine step(m,n,x,y,d,dr)
      integer x(m),y(m)
      real*4 d(m)
      dm=0.
      do 1 i=1,n
      d4=9.d9
      d3=9.d9
      d2=9.d9
      d1=9.d9
      do 2 j=1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      dd= xx*xx+yy*yy
      if(dd.gt.d4.or.dd.le.1.) go to 2
      if(dd.lt.d1) d1=dd
      if(dd.lt.d2.and.dd.gt.d1) d2=dd
      if(dd.lt.d3.and.dd.gt.d2) d3=dd
      if(dd.lt.d4.and.dd.gt.d3) d4=dd
    2 continue
      d1=sqrt(d1)
      d2=sqrt(d2)
      d3=sqrt(d3)
      d4=sqrt(d4)
    1 d(i)=(d1+d2+d3+d4)/4.0
      call dmedian(m,n,d,dr)
      return
      end
c***********************************************************************
c          Cell design 1                                               *
c***********************************************************************
      subroutine celdascapa(ax,bx,ay,by,zs,tz,avm,ms,ns,x,y,z,
     -mc,nc,xb,yb,zt,zb,dx,dy,at)
      integer xb(mc),yb(mc),zt(mc),zb(mc),dx(mc),dy(mc),
     -x(ms),y(ms),z(ms) 
      character*50 texto
      real*4 at(mc)
      zi=zs-tz                       
      zz=zs-tz/2
      tx=tz     
      ty1=tz/3  
      ty2=tz*3  
      tt=tz*tz 
      do 7 i=1,999                
      xx=ax+tx*(i-1)
      if(xx.gt.bx) go to 9
      ty=tz
      yr=ay
      do 6 j=1,999                 
      if(yr.ge.by) go to 7
      c=ty
      l=0
    3 yy=yr+c/2
      call at2(ms,ns,x,y,z,xx,yy,zz,tx,c,tz,se,am,zp)
      if(zp.lt.zs.or.am.le.0) go to 6      
      av=am/avm
      if(abs(av-1.).le.0.1) go to 4        
      c=c/av**0.41
      if(c.lt.ty1.or.c.gt.ty2) go to 6    
      l=l+1
      if(l.le.10) go to 3
      go to 6
    4 if(c.gt.9999) go to 6
      ty=c                           
      nc=nc+1
      if(nc.ge.mc) then 
        write(*,*) ' >> Warning: number of cells >',mc 
        stop   
      endif  
      xb(nc)=xx
      yb(nc)=yy
      dx(nc)=tx
      dy(nc)=ty
      zt(nc)=zs
      zb(nc)=zi
      at(nc)=av
    6 yr=yr+ty
    7 continue
    9 return
      end
c***********************************************************************
c               Cells design 2                                         * 
c***********************************************************************
      subroutine especapa(zs,avm,ms,ns,x,y,z,tz)
      integer x(ms),y(ms),z(ms)
      l=0
      c=tz                            
    1 zz=zs-c*0.5                     
      k=0
      am=-9.d9
      do 2 i=1,ns
      call at2(ms,ns,x,y,z,x(i)*1.5,y(i)*1.5,zz,c,c,c,se,aa,zp)
      if((zp-zz).le.100) go to 2
      if(aa.gt.am) am=aa
      k=k+1
    2 continue
      if(k.le.1) go to 9
      am=am/avm
      c=c/am**0.13
      l=l+1
      if(abs(am-1.).gt.0.1.and.l.lt.12) go to 1
      if(c.gt.tz) tz=c
    9 return
      end
c************************************************************************
c*************************End code **************************************