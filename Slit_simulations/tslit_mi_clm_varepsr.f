      implicit double precision (a-h,o-z)      
      include 't.inc'    
      pi = dacos(-1.d0)
      fourpi = 4.d0*pi
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      elch = 1.602D-19
      avno = 6.02214D23
      bk = 1.38066D-23
      faraday = 8.85418782D-12
      dielc = 78.3d0
      temp = 298.d0
      hpi = 0.5d0*pi
      mmm= 47
      lll= 49
      kkk = 45
      isl = 30
      open (mmm,file='input',form='formatted')
      open (kkk,file='coordm',form='formatted')
      open (lll,file='coords',form='formatted')      
      open (isl,file='slump',form='formatted')
      inzs = 69
      inzm = 76
      open(inzm, file='nzm',form='formatted')
      open(inzs, file='nzs',form='formatted')      
      rewind isl
      read(isl,*) islu
      rewind mmm
      read(mmm,*) 
      read(mmm,*) Nsub
      read(mmm,*) 
      read(mmm,*) Nm,Ns
      read(mmm,*) 
      read(mmm,*) ed2,h,arperch
      read(mmm,*) 
      read(mmm,*) dielc,cdielc
      write(*,*) 'dielc,cdielc = ',dielc,cdielc
      diffdielc = dielc-cdielc
      read(mmm,*) 
      read(mmm,*) dhspp,dhsnn
      dhsnn = dhspp      
      read(mmm,*) 
      read(mmm,*) fracwidtrams
      read(mmm,*) 
      read(mmm,*) kread,knsamp
      read(mmm,*) 
      read(mmm,*) deltaz      
      read(mmm,*) 
      read(mmm,*) dpsp,dpsn
      read(mmm,*) 
      read(mmm,*) fracclm,Rclmax,cldps
      write(*,*) 'OBS.! Varying cluster radius, Rclmax = ',Rclmax          
      write(*,*) 'cldps = ',cldps
      write(*,*) 'fracclm = ',fracclm            

      rnval = 1.d0
      rnsval = 1.d0
      rnval2 = rnval**2
      write(*,*) 'Nm,Ns = ',Nm,Ns            
      write(*,*) 'rnval = ',rnval
      write(*,*) 'rnsval = 1',rnsval

      rNm = dfloat(Nm)
      rNs = dfloat(Ns)      
      bjerrum = elch*elch/(fourpi*bk*temp*faraday*dielc*1.e-10)
      rrT = bjerrum
      write(*,*) 'bjerrum/AA = rrT/AA = ',bjerrum,rrT
      ed = 0.5d0*ed2
      edsq = ed*ed
      red2 = 1.d0/ed2
      red = 1.d0/ed

      dhs = dhspp      
      halfh = 0.5d0*h
      ih = nint(h/deltaz)
      ihalfh = ih/2
      rdeltaz = 1.d0/deltaz
      hdz = 0.5d0*deltaz
      hdhs = h-dhs
      
      dhspn = 0.5d0*(dhspp+dhsnn)      
      dhspp2 = dhspp**2
      dhspn2 = dhspn**2
      dhsnn2 = dhsnn**2

cc     rcut = 1.5d0*dhspp
c     changed 240417:      
c      rcut = 2.d0*dhspp
      rcut = dhspp+3.d0

      rcut2 = rcut**2
      rhspp = dhspp/2.d0
      redsq = 1.d0/edsq
      volume = ed2**3
      area = ed2**2
      valtot = rNm-rNs
      if (dabs(arperch).lt.100000.d0) then
      sdens = 1.d0/arperch
      else
      sdens = 0.d0
      endif
      rarea = 1.d0/area      
      write(*,*) 'total valency, valtot = ',valtot
      write(*,*) 'ed,shuss = ',ed,shuss
      write(*,*) 'centrally placed charges!!! '   
      write(*,*) 'dhspp = dhsnn = ',dhspp,dhsnn
      write(*,*) 'dhspn = dhspp = ',dhspn
      write(*,*) 'rhspp = ',rhspp
      write(*,*) 'rcut = dhspp+3\AA = ',rcut,rcut2
      write(*,*) 'epsr = cdielc+diffdielc*(r-dhspp)/dhspp, r.lt.rcut'  
      eratio = diffdielc/dhspp
      write(*,*) 'eratio,diffdielc = ',eratio,diffdielc
      write(*,*) 'Nm = ',Nm
      write(*,*) 'Nsub = ',Nsub
      write(*,*) 'lateral length of box (ed2): ',ed2
      write(*,*) 'surface separation (h): ',h
      write(*,*) 'surface charge density, sdens = ',sdens

      rNsurf = 2.d0*sdens*area
      write(*,*) 'rNsurf = ',rNsurf
c     checking electroneutrality:
      elec = rNm-rNs+rNsurf
      write(*,*) 'elec:',elec
      if (dabs(elec).gt.0.0001d0) goto 9999
      pwfact = bjerrum*sdens
      chpw = bjerrum*sdens
      wfact = -2.d0*pi*pwfact
      
      write(*,*) 'displacement parameter, dpsp = ',dpsp
      write(*,*) 'displacement parameter, dpsn = ',dpsn      
      write(*,*) 'knsamp (rdf distr.) = ',knsamp
      write(*,*) 'deltaz = ',deltaz
c      write(*,*) 'fracwid (rdf distr.) = ',fracwid      
      etam = pi*dhspp**3/6.d0
      etas = pi*dhsnn**3/6.d0      
      eta = (etam+etas)*rNm/ed2**3
      write(*,*) 'volume fraction (linear) = ',eta      
      if (kread.eq.1) then
      rewind kkk
      rewind lll
      do i = 1,Nm
      read(kkk,*) xhp(i),yhp(i),zhp(i)
      enddo
      do i = 1,Ns      
      read(lll,*) xhn(i),yhn(i),zhn(i)
      enddo

      ih2 = ih/2
      do i = 1,ih
      cz(i) = 0.d0
      ext(i) = 0.d0
      enddo
      ext(0) = 0.d0
      ext(ih+1) = 0.d0      
      write(*,*) 'reading external field'  
      rewind inzm
      rewind inzs     
      do i = 1,ih
      read(inzm,*) z,rhozm(i)
      read(inzs,*) z,rhozs(i)
      cz(i) = rhozm(i)-rhozs(i)
      enddo     
      j = ih2+1
      do i = ih2+1,ih 
      j = j-1
      cz(i) = 0.5d0*(cz(i)+cz(j))
      cz(j) = cz(i)
      enddo
      tzS6 = -1.5d0*deltaz
      do j = 0,ih+1
      tzS6 = tzS6+deltaz
      tuw = 0.d0
      ttz = -0.5d0*deltaz
      do i = 1,ih
      ttz = ttz+deltaz
      diffz = dabs(ttz-tzS6)
      tuw = cz(i)*(-2.d0*pi*diffz-Phiw(diffz))+tuw
      enddo
      ext(j) = bjerrum*tuw*deltaz
      enddo
      
      elseif (kread.eq.2) then
      do i = 1,Nm
      read(kkk,*) xhp(i),yhp(i),zhp(i)
c      read(kkk,*) t1,t2,t3
      enddo
      do i = 1,Ns      
      read(lll,*) xhn(i),yhn(i),zhn(i)
c      read(lll,*) t1,t2,t3      
      enddo
      else
      CALL cinit
      endif
      Nmac = 10
      rrNm = 1./rNm
c     convenient re-definition:
      A4 = A4/bjerrum      
      uext = 0.d0
      utot = 0.d0      
      CALL ucalc
      upart = utot
      write(*,*) 'upart = ',upart
      write(*,*) 'uext = ',uext
c      write(*,*) 'uchw = ',uchw

c      rewind inzm
c      rewind inzs
      asampnz = 0.
      do j = 1,maxnz
      rhozm(i) = 0.d0
      rhozs(i) = 0.d0
      enddo
c      ia =  90
c      open (ia,file='a.pqr',form='formatted')
c      rewind ia
c      do i = 1,Nm
c      write(ia,'(A,i7,A,i6,1f12.2,2f9.2,2f6.2)')
c     *'ATOM',i,' unk unk ',i,xhp(i),yhp(i),zhp(i),0.d0,2.d0
c      enddo
c
c     stop
      dr = deltaz
      rewind 78
      r = dhspp-dr
      nnn = nint(10./dr)+1
      do i = 1,nnn
      r = r+dr
      rmm2 = r**2
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc
      else
      tbjerrum = bjerrum
      endif               
      pot4 = tbjerrum/r
      write(78,*) r,pot4
      enddo      
c      stop

      rhphn = 0.
      attdp = 0.d0
      accdp = 0.d0
      rejdp = 0.d0
      attdn = 0.d0
      accdn = 0.d0
      rejdn = 0.d0
      attcrep = 0.
      acccrep = 0.
      rejcrep = 0.
      attdesp = 0.
      accdesp = 0.
      rejdesp = 0.
      attcren = 0.
      acccren = 0.
      rejcren = 0.
      attdesn = 0.
      accdesn = 0.
      rejdesn = 0.
      sampN = 0.
      avNm = 0.d0
      avNs = 0.d0
      avvaltot = 0.d0
      attclp = 0.
      attcln = 0.
      accclp = 0.
      acccln = 0.
      rejclp = 0.
      rejcln = 0.
      maxNn = Ns
      maxNp = Nm
      minNn = Ns
      minNp = Nm
      rdV = rarea*rdeltaz      
      do kmac = 1,Nmac
      do ksub = 1,Nsub
         
c     attempted displacement
      if (ran2(islu).lt.0.5d0) then
c     attempted displacement of cation
      ks = int(1.d0+rNm*ran2(islu))
      if (ks.gt.Nm) ks = Nm
      if (ks.lt.1) ks = 1                  
      if (ran2(islu).lt.fracclm) then
c     attempted cluster displacement of cation
      Rclust = dhspp+ran2(islu)*(Rclmax-dhspp)
      Rclust2 = Rclust*Rclust
      attclp = attclp+1.
      call clp
      if (idonk.eq.1) goto 927
      dut = (udn-udo)+def 
      if (dut.gt.100.d0) goto 927      
      if (dut.lt.0.d0) goto 9162
      if (dexp(-dut).lt.ran2(islu)) goto 927
 9162 accclp = accclp+1.      
      upart = upart+(udn-udo)
      uext = uext+def
      xhp(ks) = txs
      yhp(ks) = tys
      zhp(ks) = tzs
      do kk = 1,kcln
      km = iclustn(kk)
      xhn(km) = xcn(km)
      yhn(km) = ycn(km)
      zhn(km) = zcn(km)
      enddo
      do kk = 1,kclp
      km = iclustp(kk)
      xhp(km) = xcp(km)
      yhp(km) = ycp(km)
      zhp(km) = zcp(km)
      enddo
      goto 6
 927  rejclp = rejclp+1.
      goto 6
      else
c     attempted ordinary displacement of cation         
      attdp = attdp+1.
      call dispp
      if (idonk.eq.1) goto 203
      delta = (un-uo)+def
      if (delta.lt.0.) goto 527
      if (dexp(-delta).lt.ran2(islu)) go to 203
 527  accdp = accdp+1.
      xhp(ks) = tx
      yhp(ks) = ty
      zhp(ks) = tz
      uext = uext+def
      upart = upart+(un-uo)
      goto 6
 203  rejdp = rejdp+1.
      goto 6
      endif
      else
c     attempted displacement of anion            
      ks = int(1.d0+rNs*ran2(islu))
      if (ks.gt.Ns) ks = Ns
      if (ks.lt.1) ks = 1         

      if (ran2(islu).lt.fracclm) then
c     attempted cluster displacement of anion
      Rclust = dhspp+ran2(islu)*(Rclmax-dhspp)
      Rclust2 = Rclust*Rclust
      attcln = attcln+1.         
      call cln

      if (idonk.eq.1) goto 727
      dut = (udn-udo)+def 
      if (dut.gt.100.d0) goto 727      
      if (dut.lt.0.d0) goto 9762
      if (dexp(-dut).lt.ran2(islu)) goto 727
 9762 acccln = acccln+1.
      upart = upart+(udn-udo)
      uext = uext+def      
      xhn(ks) = txs
      yhn(ks) = tys
      zhn(ks) = tzs
      do kk = 1,kcln
      km = iclustn(kk)
      xhn(km) = xcn(km)
      yhn(km) = ycn(km)
      zhn(km) = zcn(km)
      enddo
      do kk = 1,kclp
      km = iclustp(kk)
      xhp(km) = xcp(km)
      yhp(km) = ycp(km)
      zhp(km) = zcp(km)
      enddo
      goto 6
 727  rejcln = rejcln+1.
      goto 6
      else
c     attempted ordinary displacement of anion                  
      attdn = attdn+1.
      call dispn
      if (idonk.eq.1) goto 403
      delta = (un-uo)+def 
      if (delta.lt.0.) goto 327
      if (dexp(-delta).lt.ran2(islu)) go to 403
 327  accdn = accdn+1.
      xhn(ks) = tx
      yhn(ks) = ty
      zhn(ks) = tz
      uext = uext+def
      upart = upart+(un-uo)
      goto 6
 403  rejdn = rejdn+1.
      goto 6
      endif
      endif

 6    continue
      sampN = sampN+1.
      avNs = rNs+avNs
      avNm = rNm+avNm
      avvaltot = valtot+avvaltot
      if (Ns.gt.maxNn) maxNn = Ns
      if (Nm.gt.maxNp) maxNp = Nm
      if (Ns.lt.minNn) minNn = Ns
      if (Nm.lt.minNp) minNp = Nm            
      rNmdist(Nm) = rNmdist(Nm)+1.
      rNsdist(Ns) = rNsdist(Ns)+1.            

      if (mod(ksub,knsamp).eq.0) then 
      asampnz = asampnz+1.d0
      do i = 1,Nm
      idw = int(zhp(i)*rdeltaz)+1
      rhozm(idw) = rhozm(idw)+rdV
      enddo
      do i = 1,Ns
      idw = int(zhn(i)*rdeltaz)+1      
      rhozs(idw) = rhozs(idw)+rdV
      enddo
      endif
         
      enddo
      write(*,*) 
      write(*,*) kmac
      supart = upart
      suext = uext
      uext = 0.d0
      utot = 0.d0
      CALL ucalc
      upart = utot
      dup = upart-supart
      duext = uext-suext
      write(*,*) 'upart,supart = ',upart,supart
      write(*,*) 'dup,dup/upart = ',dup,dup/upart
      write(*,*) 'uext,suext = ',uext,suext
      write(*,*) 'duext,duext/uext = ',duext,duext/uext      
      write(*,*) 'faccdp,faccdn = ',accdp/attdp,accdn/attdn
      if (attclp.gt.0.1) then
      write(*,*) 'faccclp = ',accclp/attclp
      write(*,*) attclp,accclp,rejclp
      endif
      if (attcln.gt.0.1) then
      write(*,*) 'facccln = ',acccln/attcln
      write(*,*) attcln,acccln,rejcln
      endif      
      enddo
      rewind kkk
      rewind lll
      do i = 1,Nm
      write(kkk,*) xhp(i),yhp(i),zhp(i)
      enddo
      do i = 1,Ns
      write(lll,*) xhn(i),yhn(i),zhn(i)
      enddo

      rewind inzm
      rewind inzs      
      do i = 1,ih
      densm = rhozm(i)/(asampnz)
      denss = rhozs(i)/(asampnz)
      avnm = densm+avnm
      avns = denss+avns
      apps = apps+(densm-denss)*deltaz
      if (denss.lt.1E-14) denss = 0.d0
      write(inzs,*) dfloat(i)*deltaz-0.5d0*deltaz,denss
      if (densm.lt.1E-14) densm = 0.d0
      write(inzm,*) dfloat(i)*deltaz-0.5d0*deltaz,densm
      enddo
      avnm = avnm*deltaz/h
      write(*,*) 'avnm (across h)= ',avnm
      avns = avns*deltaz/h
      write(*,*) 'avns (across h)= ',avns
      rewind isl
      write(isl,*) islu
 9999 continue
      STOP                                                              
      END


      subroutine cln
      implicit double precision (a-h,o-z)      
      include 't.inc'                              
      idonk = 0
      tef = 0.d0
      ttef = 0.d0
      dpx = cldps*(ran2(islu)-0.5d0)
      dpy = cldps*(ran2(islu)-0.5d0)
      dpz = cldps*(ran2(islu)-0.5d0)
      tzs = zhn(ks)+dpz
      if (tzs.gt.hdhs) goto 927
      if (tzs.lt.dhs) goto 927
      txs = xhn(ks)+dpx
      tys = yhn(ks)+dpy      
      if (txs.gt.ed) then
      txs = txs-ed2
      elseif (txs.lt.-ed) then
      txs = txs+ed2
      endif
      if (tys.gt.ed) then
      tys = tys-ed2
      elseif (tys.lt.-ed) then
      tys = tys+ed2
      endif
      kclp = 0
      ttxs = xhn(ks)
      ttys = yhn(ks)
      ttzs = zhn(ks)

      illt = int((tzs+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      ttef = -((zht-tzs)*ext(illt)+(tzs-zlt)*ext(illt+1))*rdeltaz
      illt = int((ttzs+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      tef = -((zht-ttzs)*ext(illt)+(ttzs-zlt)*ext(illt+1))*rdeltaz
      
      udo = 0.d0
      do k = 1,Nm
      xcp(k) = xhp(k)
      ycp(k) = yhp(k)
      zcp(k) = zhp(k)
      dx = dabs(ttxs-xhp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xhp(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = yhp(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zhp(k)+dpz
c      if (zz.gt.ed) then
c      zz = zz-ed2
c      elseif (zz.lt.-ed) then
c      zz = zz+ed2
c     endif
      if (zz.gt.hdhs) goto 927
      if (zz.lt.dhs) goto 927      
      zcp(k) = zz
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif         
      udo = -tbjerrum/dsqrt(rmm2)+udo
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)
c      udo = A4*(((r-dhspp)/dhspp)**4-1.d0)+udo
c      endif
      endif
      enddo
      
      kcln = 0
      do k = 1,ks-1
      xcn(k) = xhn(k)
      ycn(k) = yhn(k)
      zcn(k) = zhn(k)
      dx = dabs(ttxs-xhn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xhn(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = yhn(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zhn(k)+dpz
c      if (zz.gt.ed) then
c      zz = zz-ed2
c      elseif (zz.lt.-ed) then
c      zz = zz+ed2
c     endif
      if (zz.gt.hdhs) goto 927
      if (zz.lt.dhs) goto 927            
      zcn(k) = zz
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udo = tbjerrum/dsqrt(rmm2)+udo
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udo = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udo
c      endif
      endif
      enddo

      do k = ks+1,Ns
      xcn(k) = xhn(k)
      ycn(k) = yhn(k)
      zcn(k) = zhn(k)
      dx = dabs(ttxs-xhn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xhn(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = yhn(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zhn(k)+dpz
c      if (zz.gt.ed) then
c      zz = zz-ed2
c      elseif (zz.lt.-ed) then
c      zz = zz+ed2
c     endif
      if (zz.gt.hdhs) goto 927
      if (zz.lt.dhs) goto 927            
      zcn(k) = zz
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udo = tbjerrum/dsqrt(rmm2)+udo
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udo = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udo
c      endif         
      endif
      enddo

      udn = 0.d0
      kkclp = 0
      do k = 1,Nm
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udn = -tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif         
      endif   
      enddo

      kkcln = 0
      do k = 1,ks-1
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif         
      endif   
      enddo
      do k = ks+1,Ns
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif         
      endif
      enddo
      
      if (kkcln.ne.kcln.or.kkclp.ne.kclp) goto 927      
c     no cluster overlaps => calculate 
c     present and trial energies with other cluster members
      xcn(ks) = txs
      ycn(ks) = tys
      zcn(ks) = tzs      
      do kk = 1,kclp
      km = iclustp(kk)
      stx = xcp(km)
      sty = ycp(km)
      stz = zcp(km)
      ttx = xhp(km)
      tty = yhp(km)
      ttz = zhp(km)

      illt = int((stz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      ttef = ((zht-stz)*ext(illt)+(stz-zlt)*ext(illt+1))*rdeltaz+ttef
      illt = int((ttz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      tef = ((zht-ttz)*ext(illt)+(ttz-zlt)*ext(illt+1))*rdeltaz+tef     
      
      do k = 1,km-1
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = tbjerrum/dsqrt(rmm2)+udn      
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upp(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo      
      enddo
      do k = km+1,Nm
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upp(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo            
      enddo

      do k = 1,Ns      
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = -tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upn(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo       
      enddo
      enddo

      do kk = 1,kcln
      km = iclustn(kk)
      stx = xcn(km)
      sty = ycn(km)
      stz = zcn(km)
      ttx = xhn(km)
      tty = yhn(km)
      ttz = zhn(km)

      illt = int((stz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      ttef = -((zht-stz)*ext(illt)+(stz-zlt)*ext(illt+1))*rdeltaz+ttef
      illt = int((ttz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      tef = -((zht-ttz)*ext(illt)+(ttz-zlt)*ext(illt+1))*rdeltaz+tef     
      
      do k = 1,km-1
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                     
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upp(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo      
      enddo
      do k = km+1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upp(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo            
      enddo

      do k = 1,Nm
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = -tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upn(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo       
      enddo
      enddo
      def = ttef-tef      
      goto 929
 927  idonk = 1
 929  continue
      return
      end            

      subroutine clp
      implicit double precision (a-h,o-z)      
      include 't.inc'                              
      idonk = 0
      dpx = cldps*(ran2(islu)-0.5d0)
      dpy = cldps*(ran2(islu)-0.5d0)
      dpz = cldps*(ran2(islu)-0.5d0)
      tzs = zhp(ks)+dpz
      if (tzs.gt.hdhs) goto 927
      if (tzs.lt.dhs) goto 927
      txs = xhp(ks)+dpx
      tys = yhp(ks)+dpy
      if (txs.gt.ed) then
      txs = txs-ed2
      elseif (txs.lt.-ed) then
      txs = txs+ed2
      endif
      if (tys.gt.ed) then
      tys = tys-ed2
      elseif (tys.lt.-ed) then
      tys = tys+ed2
      endif

      kcln = 0
      ttxs = xhp(ks)
      ttys = yhp(ks)
      ttzs = zhp(ks)

      illt = int((tzs+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      ttef = ((zht-tzs)*ext(illt)+(tzs-zlt)*ext(illt+1))*rdeltaz
      illt = int((ttzs+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      tef = ((zht-ttzs)*ext(illt)+(ttzs-zlt)*ext(illt+1))*rdeltaz
      
      udo = 0.d0
      do k = 1,Ns
      xcn(k) = xhn(k)
      ycn(k) = yhn(k)
      zcn(k) = zhn(k)
      dx = dabs(ttxs-xhn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xhn(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = yhn(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zhn(k)+dpz
c      if (zz.gt.ed) then
c      zz = zz-ed2
c      elseif (zz.lt.-ed) then
c      zz = zz+ed2
c     endif
      if (zz.gt.hdhs) goto 927
      if (zz.lt.dhs) goto 927            
      zcn(k) = zz
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udo = -tbjerrum/dsqrt(rmm2)+udo
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udo = A4*(((r-dhspp)/dhspp)**4-1.d0)+udo
c      endif         
      endif
      enddo
      
      kclp = 0
      do k = 1,ks-1
      xcp(k) = xhp(k)
      ycp(k) = yhp(k)
      zcp(k) = zhp(k)
      dx = dabs(ttxs-xhp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xhp(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = yhp(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zhp(k)+dpz
c      if (zz.gt.ed) then
c      zz = zz-ed2
c      elseif (zz.lt.-ed) then
c      zz = zz+ed2
c     endif
      if (zz.gt.hdhs) goto 927
      if (zz.lt.dhs) goto 927            
      zcp(k) = zz
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udo = tbjerrum/dsqrt(rmm2)+udo
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udo = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udo
c      endif         
      endif
      enddo

      do k = ks+1,Nm
      xcp(k) = xhp(k)
      ycp(k) = yhp(k)
      zcp(k) = zhp(k)
      dx = dabs(ttxs-xhp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xhp(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = yhp(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zhp(k)+dpz
c      if (zz.gt.ed) then
c      zz = zz-ed2
c      elseif (zz.lt.-ed) then
c      zz = zz+ed2
c     endif
      if (zz.gt.hdhs) goto 927
      if (zz.lt.dhs) goto 927            
      zcp(k) = zz
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udo = tbjerrum/dsqrt(rmm2)+udo
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udo = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udo
c      endif         
      endif
      enddo            

      udn = 0.d0
      kkcln = 0
      do k = 1,Ns
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udn = -tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif         
      endif   
      enddo

      kkclp = 0
      do k = 1,ks-1
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif         
      endif   
      enddo
      do k = ks+1,Nm
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif                  
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif         
      endif   
      enddo
      
      if (kkcln.ne.kcln.or.kkclp.ne.kclp) goto 927      
c     no cluster overlaps => calculate 
c     present and trial energies with other cluster members
      xcp(ks) = txs
      ycp(ks) = tys
      zcp(ks) = tzs      
      do kk = 1,kcln
      km = iclustn(kk)
      stx = xcn(km)
      sty = ycn(km)
      stz = zcn(km)
      ttx = xhn(km)
      tty = yhn(km)
      ttz = zhn(km)

      illt = int((stz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      ttef = -((zht-stz)*ext(illt)+(stz-zlt)*ext(illt+1))*rdeltaz+ttef
      illt = int((ttz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      tef = -((zht-ttz)*ext(illt)+(ttz-zlt)*ext(illt+1))*rdeltaz+tef           
      
      do k = 1,km-1
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upp(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo      
      enddo
      do k = km+1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upp(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo            
      enddo

      do k = 1,Nm      
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = -tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upn(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo       
      enddo
      enddo

      do kk = 1,kclp
      km = iclustp(kk)
      stx = xcp(km)
      sty = ycp(km)
      stz = zcp(km)
      ttx = xhp(km)
      tty = yhp(km)
      ttz = zhp(km)

      illt = int((stz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      ttef = ((zht-stz)*ext(illt)+(stz-zlt)*ext(illt+1))*rdeltaz+ttef
      illt = int((ttz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      tef = ((zht-ttz)*ext(illt)+(ttz-zlt)*ext(illt+1))*rdeltaz+tef     
      
      do k = 1,km-1
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upp(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo      
      enddo
      do k = km+1,Nm
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = +A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upp(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo            
      enddo

      do k = 1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      udn = -tbjerrum/dsqrt(rmm2)+udn
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      udn = A4*(((r-dhspp)/dhspp)**4-1.d0)+udn
c      endif      
      udo = upn(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo       
      enddo
      enddo

      def = ttef-tef      
      
      goto 929
 927  idonk = 1
 929  continue
      return
      end      
      
      subroutine cinit
      implicit double precision (a-h,o-z)      
      include 't.inc'                        
      do i = 1,Nm
      xhp(i) = (ran2(islu)-0.5d0)*ed2
      yhp(i) = (ran2(islu)-0.5d0)*ed2
      zhp(i) = dhs+ran2(islu)*(h-2.d0*dhs)            
      xhn(i) = (ran2(islu)-0.5d0)*ed2
      yhn(i) = (ran2(islu)-0.5d0)*ed2
      zhn(i) = dhs+ran2(islu)*(h-2.d0*dhs)                
      enddo
      write(*,*) 'Ions successfully inserted'                            
      return
      end

      subroutine dispp
      implicit double precision (a-h,o-z)
      include 't.inc' 
      idonk = 0
      if (ran2(islu).lt.0.1d0) then
      tz = dhs+ran2(islu)*(h-2.d0*dhs) 
      tx = ed2*(ran2(islu)-0.5d0)
      ty = ed2*(ran2(islu)-0.5d0)
      else         
      dtz = dpsp*(ran2(islu)-0.5d0)
      tz = zhp(ks)+dtz
      if (tz.gt.hdhs) goto 127
      if (tz.lt.dhs) goto 127         
      dtx = dpsp*(ran2(islu)-0.5d0)
      dty = dpsp*(ran2(islu)-0.5d0)
      tx = xhp(ks)+dtx
      if (tx.gt.ed) then
      tx = tx-ed2
      elseif (tx.lt.-ed) then
      tx = tx+ed2
      endif
      ty = yhp(ks)+dty
      if (ty.gt.ed) then
      ty = ty-ed2
      elseif (ty.lt.-ed) then
      ty = ty+ed2
      endif
      endif
      
c     do k = 1,Nm
      do k = 1,Ns      
      idonk = kdonk(tx,ty,tz,xhn(k),yhn(k),zhn(k),dhspn2)
      if (idonk.eq.1) goto 98
      enddo
      do k = 1,ks-1
      idonk = kdonk(tx,ty,tz,xhp(k),yhp(k),zhp(k),dhspp2)
      if (idonk.eq.1) goto 98
      enddo
      do k = ks+1,Nm
      idonk = kdonk(tx,ty,tz,xhp(k),yhp(k),zhp(k),dhspp2)
      if (idonk.eq.1) goto 98
      enddo

      un = 0.d0
      uo = 0.d0
      qttx = xhp(ks)
      qtty = yhp(ks)
      qttz = zhp(ks)
      do k = 1,ks-1
      un = upp(tx,ty,tz,xhp(k),yhp(k),zhp(k))+un
      uo = upp(qttx,qtty,qttz,xhp(k),yhp(k),zhp(k))+uo      
      enddo
      do k = ks+1,Nm
      un = upp(tx,ty,tz,xhp(k),yhp(k),zhp(k))+un
      uo = upp(qttx,qtty,qttz,xhp(k),yhp(k),zhp(k))+uo      
      enddo
      do k = 1,Ns
      un = upn(tx,ty,tz,xhn(k),yhn(k),zhn(k))+un
      uo = upn(qttx,qtty,qttz,xhn(k),yhn(k),zhn(k))+uo      
      enddo

      illt = int((tz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      ttef = ((zht-tz)*ext(illt)+(tz-zlt)*ext(illt+1))*rdeltaz
c      tuch = wfact*tz+wfact*(h-tz)
      illt = int((qttz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      tef = ((zht-qttz)*ext(illt)+(qttz-zlt)*ext(illt+1))*rdeltaz
c      uch = wfact*qttz+wfact*(h-qttz)
c      duchw = tuch-uch      
      def = ttef-tef
      
      goto 98
 127  idonk = 1
 98   continue
      return
      end

      subroutine dispn
      implicit double precision (a-h,o-z)
      include 't.inc' 
      idonk = 0
      if (ran2(islu).lt.0.1d0) then
      tz = dhs+ran2(islu)*(h-2.d0*dhs)          
      tx = ed2*(ran2(islu)-0.5d0)
      ty = ed2*(ran2(islu)-0.5d0)
      else               
      dtz = dpsn*(ran2(islu)-0.5d0)
      tz = zhn(ks)+dtz         
      if (tz.gt.hdhs) goto 127
      if (tz.lt.dhs) goto 127                  
      dtx = dpsn*(ran2(islu)-0.5d0)
      dty = dpsn*(ran2(islu)-0.5d0)      
      tx = xhn(ks)+dtx
      if (tx.gt.ed) then
      tx = tx-ed2
      elseif (tx.lt.-ed) then
      tx = tx+ed2
      endif
      ty = yhn(ks)+dty
      if (ty.gt.ed) then
      ty = ty-ed2
      elseif (ty.lt.-ed) then
      ty = ty+ed2
      endif
      endif

      do k = 1,Nm
      idonk = kdonk(tx,ty,tz,xhp(k),yhp(k),zhp(k),dhspn2)
      if (idonk.eq.1) goto 98
      enddo
      do k = 1,ks-1
      idonk = kdonk(tx,ty,tz,xhn(k),yhn(k),zhn(k),dhsnn2)
      if (idonk.eq.1) goto 98
      enddo
      do k = ks+1,Ns      
      idonk = kdonk(tx,ty,tz,xhn(k),yhn(k),zhn(k),dhsnn2)
      if (idonk.eq.1) goto 98
      enddo

      un = 0.d0
      uo = 0.d0
      qttx = xhn(ks)
      qtty = yhn(ks)
      qttz = zhn(ks)
      do k = 1,Nm
      un = upn(tx,ty,tz,xhp(k),yhp(k),zhp(k))+un
      uo = upn(qttx,qtty,qttz,xhp(k),yhp(k),zhp(k))+uo      
      enddo
      do k = 1,ks-1
      un = upp(tx,ty,tz,xhn(k),yhn(k),zhn(k))+un
      uo = upp(qttx,qtty,qttz,xhn(k),yhn(k),zhn(k))+uo      
      enddo
      do k = ks+1,Ns      
      un = upp(tx,ty,tz,xhn(k),yhn(k),zhn(k))+un
      uo = upp(qttx,qtty,qttz,xhn(k),yhn(k),zhn(k))+uo      
      enddo

      illt = int((tz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      ttef = -((zht-tz)*ext(illt)+(tz-zlt)*ext(illt+1))*rdeltaz
c      tuch = -wfact*tz-wfact*(h-tz)
      illt = int((qttz+hdz)*rdeltaz)
      zlt = dfloat(illt)*deltaz-hdz
      zht = zlt+deltaz
      tef = -((zht-qttz)*ext(illt)+(qttz-zlt)*ext(illt+1))*rdeltaz
c      uch = -wfact*qttz-wfact*(h-qttz)
c      duchw = tuch-uch            
      def = ttef-tef      
      goto 98
 127  idonk = 1
      
 98   continue
      return
      end

      subroutine ucalc
      implicit double precision (a-h,o-z)
      include 't.inc'
      uext = 0.d0
      utot = 0.d0 
      do i = 1,Nm
      tx = xhp(i)
      ty = yhp(i)
      tz = zhp(i)
      if (abs(tx).gt.ed) write(*,*) '|XM|.GT.ED!',i,tx
      if (abs(ty).gt.ed) write(*,*) '|YM|.GT.ED!',i,ty
      if (tz.gt.hdhs) write(*,*)'ZM.GT.HDHS !!! ',i,tz
      if (tz.lt.dhs) write(*,*) 'ZM.LT.DHS !!! ',i,tz
      do k = i+1,Nm
      dx = abs(tx-xhp(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhp(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhp(k))
c      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) write(*,*) 'RMM2.LT.DHSPP2!',i,k,rmm2
      enddo
      do k = 1,Ns      
      dx = abs(tx-xhn(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhn(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhn(k))
c      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) write(*,*) 'RMM2.LT.DHSPN2!',i,k,rmm2
      enddo
      enddo

      do i = 1,Nm
      tx = xhp(i)
      ty = yhp(i)
      tz = zhp(i)
      do k = i+1,Nm
      dx = abs(tx-xhp(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhp(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhp(k))
c      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      r = dsqrt(rmm2)
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      utot = tbjerrum/r+utot
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      utot = +A4*(((r-dhspp)/dhspp)**4-1.d0)+utot
c      endif      
      enddo      
      do k = 1,Ns      
      dx = abs(tx-xhn(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhn(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhn(k))
c      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      r = dsqrt(rmm2)
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      utot = -tbjerrum/r+utot
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      utot = A4*(((r-dhspp)/dhspp)**4-1.d0)+utot
c      endif      
      enddo
      enddo

      do i = 1,Ns      
      tx = xhn(i)
      ty = yhn(i)
      tz = zhn(i)
      if (abs(tx).gt.ed) write(*,*) '|XS|.GT.ED!',i,tx
      if (abs(ty).gt.ed) write(*,*) '|YS|.GT.ED!',i,ty
      if (tz.gt.hdhs) write(*,*)'ZS.GT.HDHS !!! ',i,tz
      if (tz.lt.dhs) write(*,*) 'ZS.LT.DHS !!! ',i,tz
      do k = i+1,Ns      
      dx = abs(tx-xhn(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhn(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhn(k))
c      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) write(*,*) 'RMM2.LT.DHSNN2!',i,k,rmm2
      enddo
      enddo

c      stop
      
      do i = 1,Ns      
      tx = xhn(i)
      ty = yhn(i)
      tz = zhn(i)
       do k = i+1,Ns      
      dx = abs(tx-xhn(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhn(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhn(k))
c      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      r = dsqrt(rmm2)
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      utot = tbjerrum/r+utot
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      utot = +A4*(((r-dhspp)/dhspp)**4-1.d0)+utot
c      endif      
      enddo
      enddo

      uext = 0.d0
c      write(*,*) 'hdz,rdeltaz = ',hdz,rdeltaz
      do i = 1,Nm
      tzS6 = zhp(i)
      ill = int((tzS6+hdz)*rdeltaz)
      zl = dfloat(ill)*deltaz-hdz
      zh = zl+deltaz
      uext = ((zh-tzS6)*ext(ill)+(tzS6-zl)*ext(ill+1))*rdeltaz+uext
      enddo
      do i = 1,Ns
      tzS6 = zhn(i)
      ill = int((tzS6+hdz)*rdeltaz)
      zl = dfloat(ill)*deltaz-hdz
      zh = zl+deltaz
      uext = -((zh-tzS6)*ext(ill)+(tzS6-zl)*ext(ill+1))*rdeltaz+uext
      enddo
c      write(*,*) 'zl,zh,ill = ',zl,zh,ill
c      write(*,*) 'ext(ill),uext = ',ext(ill),uext
c      stop
      
      return
      end

      function kdonk(xx,yy,zz,xxx,yyy,zzz,dref2)
      implicit double precision (a-h,o-z)
      include 't.inc'
      dx = abs(xx-xxx)
      dx = dx-aint(dx/ed)*ed2
      dy = abs(yy-yyy)
      dy = dy-aint(dy/ed)*ed2
      dz = abs(zz-zzz)
c      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dref2) then
      kdonk = 1
      else
      kdonk = 0
      endif
      return
      end

      function upp(xx,yy,zz,xxx,yyy,zzz)
      implicit double precision (a-h,o-z)
      include 't.inc'
      dx = abs(xx-xxx)
      dx = dx-aint(dx*red)*ed2
      dy = abs(yy-yyy)
      dy = dy-aint(dy*red)*ed2
      dz = abs(zz-zzz)
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      upp = tbjerrum/dsqrt(rmm2)
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      upp = +A4*(((r-dhspp)/dhspp)**4-1.d0)+uc
c      else
c      upp = uc
c      endif   
      return
      end
      
      function upn(xx,yy,zz,xxx,yyy,zzz)
      implicit double precision (a-h,o-z)
      include 't.inc'
      dx = abs(xx-xxx)
      dx = dx-aint(dx*red)*ed2
      dy = abs(yy-yyy)
      dy = dy-aint(dy*red)*ed2
      dz = abs(zz-zzz)
c      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      tbjerrum = bjerrum*dielc/tdielc      
      else
      tbjerrum = bjerrum
      endif               
      upn = -tbjerrum/dsqrt(rmm2)
c      if (rmm2.lt.rcut2) then
c      r = dsqrt(rmm2)         
c      upn = A4*(((r-dhspp)/dhspp)**4-1.d0)+uc
c      else
c      upn = uc
c      endif            
      return
      end

      function dist(xx,yy,zz,xxx,yyy,zzz)
      implicit double precision (a-h,o-z)
      include 't.inc'
      dx = abs(xx-xxx)
      dx = dx-aint(dx/ed)*ed2
      dy = abs(yy-yyy)
      dy = dy-aint(dy/ed)*ed2
      dz = abs(zz-zzz)
c      dz = dz-aint(dz/ed)*ed2
      dist = dsqrt(dx*dx+dy*dy+dz*dz)
      return
      end

      function Phiw(z)
      implicit double precision (a-h,o-z)
      include 't.inc'      
      zsq = z*z
      Phiw = 
     *8.d0*ed*dlog((dsqrt(2.d0*edsq+zsq)+ed)/(dsqrt(edsq+zsq)))-
     *2.d0*z*(dasin((edsq**2-zsq*zsq-2.d0*edsq*zsq)/(edsq+zsq)**2)+
     *hpi)
      return
      end







